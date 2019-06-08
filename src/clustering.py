from collections import Counter
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.cluster import KMeans, Birch
from sklearn.mixture import GaussianMixture
from sklearn.manifold import TSNE
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from logzero import logger
from BITS.run import run_edlib
from BITS.utils import run_command, print_log, NoDaemonPool, save_pickle, load_pickle
from BITS.plot import generate_layout, show_plot, generate_scatter
from BITS.scheduler import Scheduler
import consed
#from .dpmm import DPMM


class Clustering:
    """
    Root class of clusterings of
      1. DNA sequences (regardless of cyclic alignment and/or revcomp)
      2. variant matrix (by Consed)

    <self.assignment> is the final result of the clustering.

    Both <self.s_dist_mat> (squre distance matrix) and <self.c_dist_mat> (condensed distance matrix)
    must be set when the distance matrix is computed.
    """

    def __init__(self, data, names=None):
        self.data = data
        self.N = data.shape[0]   # num of data. NOTE: Each row in <data> must be each datum
        if names is None:
            self.names = list(range(self.N))
        else:
            self.names = names   # name for each datum. used for plots
        self.assignment = np.full(self.N, -1, dtype='int8')   # cluster assignment for each data
        self.precomputed = {}   # used to avoid re-calculation of clustering results

    def _hierarchical(self, method):
        # Since computation of linkage is somewhat heavy, save it for the next time
        if method in self.precomputed:
            self.hc_result = self.precomputed[method]
        else:
            assert hasattr(self, "c_dist_mat"), "No condensed distance matrix"
            self.hc_result = linkage(self.c_dist_mat, method=method)
            self.precomputed[method] = self.hc_result

    def dendrogram(self, method="ward", figsize=(18, 10)):
        """
        Show dendrogram of data with <method>.
        """

        self._hierarchical(method)
        plt.figure(figsize=figsize)
        dendrogram(self.hc_result)
        plt.show()

    def cluster_hierarchical(self, method="ward", criterion="distance", threshold=0.7):
        """
        <method> = {single, complete, average, weighted, centroid, median, ward (default)}
        <criterion> = {inconsistent, distance (default), maxclust, monocrit, maxclust_monocrit}
        <threshold> = threshold in distance criterion

        <self.c_dist_mat> must be precomputed. Currently supports only distance criterion.
        """

        self._hierarchical(method)

        # Calculate cluster assignment for each data
        self.assignment = np.array(fcluster(self.hc_result, t=threshold, criterion=criterion))

    @property
    def n_clusters(self):
        return len(set(self.assignment))

    def cluster(self, cluster_id, return_where=False):
        """
        Return the data belonging to the cluster whose ID is specified.
        If <return_where>, return the indices of the data instead of the data itself.
        """

        where = np.where(self.assignment == cluster_id)[0]
        if return_where:
            return where
        else:
            return self.data[where]

    def clusters(self, return_where=False):
        """
        From <self.assignment>, generate each cluster with data belonging to it.
        If <return_where>, return the indices of the data instead of the data itself.
        """

        for cluster_id in sorted(list(set(self.assignment))):
            yield (cluster_id, self.cluster(cluster_id, return_where))

    def merge_cluster(self, cluster_to, cluster_from):
        self.assignment[self.cluster(cluster_from, return_where=True)] = cluster_to

    def hist_cluster_size(self, bins=100, n_print=5):
        """
        Histogram of the size of the clusters.
        Also both largest and smallest <n_print> clusters are shown.
        """

        cluster_size = Counter(self.assignment)
        print(f"{n_print} largest clusters:")
        for cluster_id, size in cluster_size.most_common()[:n_print]:
            print(f"cluster {cluster_id} ({size} units)")
        print(f"{n_print} smallest clusters:")
        for cluster_id, size in cluster_size.most_common()[-n_print:]:
            print(f"cluster {cluster_id} ({size} units)")

        plt.hist(list(cluster_size.values()), bins=bins)
        plt.show()

    def plot_dist_mat(self, show_scale=False, zmin=0, zmax=1, width=500, height=500, title=None, out_fname=None):
        """
        Draw a heatmap of the (squared) distance matrix.
        """

        #assert self.N <= 100, "Too many data to plot"

        trace = go.Heatmap(z=self.s_dist_mat,
                           colorscale="YlGnBu",
                           zmin=zmin,
                           zmax=zmax,
                           showscale=show_scale)
        layout = generate_layout(width, height, title=title)
        layout["yaxis"] = dict(autorange="reversed")
        show_plot([trace], layout, out_fname=out_fname)

    def plot_tsne(self, width=700, height=700, title=None, out_fname=None):
        """
        Embed data into a two dimensional space using t-SNE.
        <self.s_dist_mat> must be precomputed.
        """

        assert hasattr(self, "s_dist_mat"), "No distance matrix"

        coord = TSNE(n_components=2, metric='precomputed').fit_transform(self.s_dist_mat)
        trace = generate_scatter(coord[:, 0].T,
                                 coord[:, 1].T,
                                 text=[f"{self.names[i]}<br>{self.assignment[i]}"
                                       for i in range(self.N)],
                                 marker_size=3)
        trace["marker"].update(dict(color=self.assignment,
                                    colorscale="Rainbow",
                                    showscale=False))
        layout = generate_layout(width, height, title=title)
        show_plot([trace], layout, out_fname=out_fname)


def __calc_dist_array(i, data, cyclic, rc):
    """
    Compute row i vs columns (i+1) to N in the distance matrix. N is data size.
    """

    # row <i> in the distance matrix
    dist_array = np.array([run_edlib(data[i],
                                     data[j],
                                     "global",
                                     cyclic=cyclic,
                                     rc=rc,
                                     only_diff=True)
                           for j in range(i + 1, data.shape[0])],
                          dtype='float32')

    #logger.debug(f"Finished @ row {i}")
    return (i, dist_array)


def _calc_dist_array(rows, data, cyclic, rc):
    logger.debug(f"Starting row {rows[0]}-{rows[-1]}")
    return [__calc_dist_array(row, data, cyclic, rc) for row in rows]


class ClusteringSeqs(Clustering):
    """
    Perform clustering of DNA sequences given.
    Both greedy approach based on cluster diameter and hierarchical clustering approach are available.

    Run <self.cluster_greedy> for former, and <self.cluster_hierarchical> for latter.
    """

    def __init__(self, data, names=None, cyclic=True, rc=True):
        super().__init__(data, names)
        self.cyclic = cyclic   # do cyclic alignment between two sequences
        self.rc = rc   # allow reverse complement when taking alignment

    def _calc_dist_mat(self, rows, n_core):
        """
        Calculate all-vs-all distance matrix between the sequences.
        Both cyclic alignment and reverse complement sequence are considered.
        Computation is performed in parallel, and the unit of parallelization
        is each row of the (triangular) distance matrix.
        """

        n_sub = -(-len(rows) // n_core)
        tasks = [(rows[i * n_sub:(i + 1) * n_sub],
                  self.data,
                  self.cyclic,
                  self.rc)
                 for i in range(n_core)]

        dist_arrays = []
        with NoDaemonPool(n_core) as pool:
            for ret in pool.starmap(_calc_dist_array, tasks):
                dist_arrays += ret
        return dist_arrays

    @print_log("distance matrix calculation")
    def calc_dist_mat(self, n_core=1, n_distribute=1, dir_prefix=".", file_prefix="clustering"):
        """
        <prefix> must be that of only FILE NAMES, do not include directory name.
        """

        rows = [int(i / 2) if i % 2 == 0
                else self.N - 2 - int((i - 1) / 2)
                for i in range(self.N - 1)]

        if n_core * n_distribute > len(rows):
            n_core, n_distribute = 1, 1

        if n_distribute == 1:
            # No distributed computation. Directly calculate in parallel.
            dist_arrays = self._calc_dist_mat(rows, n_core)
        else:
            # Distributed computation
            n_sub = -(-len(rows) // n_distribute)
            rows_sub = [rows[i * n_sub:(i + 1) * n_sub]
                        for i in range(n_distribute)]
            n_digit = int(np.log10(n_distribute) + 1)
            jids = []
            s = Scheduler("sge",
                          "qsub",
                          out_log=f"{dir_prefix}/log")
            save_pickle(self, f"{dir_prefix}/{file_prefix}_obj.pkl")
            for i, rs in enumerate(rows_sub):
                index = str(i + 1).zfill(n_digit)
                save_pickle(rs, f"{dir_prefix}/{file_prefix}_rows.{index}.pkl")
                jids.append(s.submit(' '.join([f"calc_dist_mat_sub.py",
                                               f"-n {n_core}",
                                               f"{dir_prefix}/{file_prefix}_obj.pkl",
                                               f"{dir_prefix}/{file_prefix}_rows.{index}.pkl",
                                               f"{dir_prefix}/{file_prefix}.{index}.pkl"]),
                                     f"{dir_prefix}/{file_prefix}.sge.{index}",
                                     job_name="dist_mat_sub",
                                     n_core=n_core))

            # Merge the results
            s.submit("sleep 1s",
                     f"{dir_prefix}/gather.sge",
                     job_name="gather_calc_dist_mat",
                     depend=jids,
                     wait=True)

            dist_arrays = []
            for fname in run_command(f"find {dir_prefix} -name '{file_prefix}.*.pkl'").strip().split('\n'):
                dist_arrays += load_pickle(fname)

        self.s_dist_mat = np.zeros((self.N, self.N), dtype='float32')
        for r in dist_arrays:
            i, dist_array = r
            self.s_dist_mat[i, i + 1:] = self.s_dist_mat[i + 1:, i] = dist_array
        self.c_dist_mat = squareform(self.s_dist_mat)

    def _generate_consensus(self):
        ret = {}
        index = 0
        for cluster_id, seqs in super().clusters():
            # TODO: XXX: here assumes the first sequence is somewhat "accurate" though this does not always hold; choose the "center" sequence!
            cons_seq = consed.consensus(list(seqs) if not (self.cyclic or self.rc)   # input are already synchronized
                                        else [seq if i == 0
                                              else run_edlib(seqs.iloc[0],
                                                             seq,
                                                             "global",
                                                             cyclic=self.cyclic,
                                                             rc=self.rc,
                                                             return_seq=True).seq
                                              for i, seq in enumerate(seqs)],
                                        n_iter=3)
            if cons_seq != "":
                ret[index] = (cluster_id, seqs.shape[0], len(cons_seq), cons_seq)
                index += 1
        return pd.DataFrame.from_dict(ret,
                                      orient="index",
                                      columns=("cluster_id", "cluster_size", "length", "sequence"))

    def _merge_clusters(self, th_merge):
        flag_next = False
        n_cons = self.cons_seqs.shape[0]   # != <self.n_clusters> due to Consed error
        for i in range(n_cons - 1):
            for j in range(i + 1, n_cons):
                if run_edlib(self.cons_seqs["sequence"].iloc[i],
                             self.cons_seqs["sequence"].iloc[j],
                             "global",
                             cyclic=self.cyclic,
                             rc=self.rc,
                             only_diff=True) < th_merge:
                    self.merge_cluster(self.cons_seqs["cluster_id"].iloc[i],
                                       self.cons_seqs["cluster_id"].iloc[j])
                    flag_next = True
        self.cons_seqs = self._generate_consensus()
        return flag_next

    def generate_consensus(self,
                           th_merge=0.05,
                           th_noisy=0.01,
                           th_synchronize=0.3):
        """
        Calculate a consensus sequence for each cluster.
        1. Initial consensus sequences are computed by Consed.
        2. Every two consensus sequences which have dis-similarity less than <th_merge> are merged.
           If you do not want to merge any clusters, then specify <th_merge=0>.
        3. Any cluster whose size is less than (<self.N> * <th_noisy>) is removed.
           Specify <th_noisy=0> to keep all clusters.
        4. Synchronize the phase of clusters which share dis-similarity less than <th_synchronize>.
        """

        # TODO: detect and remove STR-like consensus sequences (in peak?)

        # Initial consensus sequences
        self.cons_seqs = self._generate_consensus()
        logger.info(f"Initial consensus sequences:\n{self.cons_seqs}")

        # Merge too close clusters
        while self._merge_clusters(th_merge):
            pass
        logger.info(f"After merging similar units:\n{self.cons_seqs}")

        # Remove remaining noisy clusters
        del_row = [index for index, df in self.cons_seqs.iterrows()
                   if df["cluster_size"] < self.N * th_noisy]   # too small cluster
        self.cons_seqs = self.cons_seqs.drop(del_row).reset_index(drop=True)
        logger.debug(f"After removing noisy clusters:\n{self.cons_seqs}")

        # Synchronize phase of the consensus units (= start position)
        n_cons = self.cons_seqs.shape[0]
        for i in range(n_cons - 1):   # TODO: simultaneously synchronize, or fix single seed
            for j in range(i + 1, n_cons):
                align = run_edlib(self.cons_seqs["sequence"].iloc[i],
                                  self.cons_seqs["sequence"].iloc[j],
                                  "global",
                                  cyclic=self.cyclic,
                                  rc=self.rc,
                                  return_seq=True)
                if align.diff < th_synchronize:
                    logger.debug(f"Synchronize {i} and {j} (strand = {align.strand})")
                    self.cons_seqs.loc[j, "sequence"] = align.seq
        logger.info(f"Final consensus sequences:\n{self.cons_seqs}")


class ClusteringNumeric(Clustering):
    """
    A child class which additionally has clustering methods only for numerical data.
    """

    def calc_dist_mat(self, metric="euclidean"):
        # NOTE: Clustering itself does not require a distance matrix, but t-SNE plot does.
        self.c_dist_mat = pdist(self.data, metric=metric)
        self.s_dist_mat = squareform(self.c_dist_mat)

    def do_clustering(self, method, **kwargs):
        """
        Wrapper of the clusterings with the data (not distance matrix).

        As <kwargs>, specify
          - "n_clusters" for "kmeans", "gmm", "birch", "nmf" (mandatory)
          - "max_cluster" for "gmm_bic" (optional)
          - "max_iteration" for "dpmm" (optional)
        """

        assert method in ("kmeans", "gmm", "birch", "nmf", "gmm_bic", "dpmm")

        params = tuple(sorted(kwargs.items()))
        if (method, params) in self.precomputed:
            self.assignment = self.precomputed[(method, params)]
        else:
            self.assignment = getattr(self, method)(**kwargs)
            self.precomputed[(method, params)] = self.assignment

    ## ------------------------------------------ ##
    ##  Clustering with fixed number of clusters  ##
    ## ------------------------------------------ ##

    def kmeans(self, n_clusters):
        return KMeans(n_clusters=n_clusters).fit_predict(self.data)

    def gmm(self, n_clusters):
        gmm = GaussianMixture(n_components=n_clusters, covariance_type='full')
        gmm.fit(self.data)
        return gmm.predict(self.data)

    def birch(self, n_clusters):
        return Birch(n_clusters=n_clusters).fit_predict(self.data)

    def nmf(self, n_clusters):   # TODO: debug
        nmf = NMF(n_components=n_clusters)
        W = nmf.fit_transform(self.data)
        H = nmf.components_
        print(W)
        print(H)

    ## ---------------------------------------------- ##
    ##  Clustering with arbitrary number of clusters  ##
    ## ---------------------------------------------- ##

    def dpmm(self, max_iteration=1000000):   # TODO: tune parameters (max_iteration, converge_threshold_*)
        """
        Dirichlet Process (Multinomial) Mixture Model with concentration hyperparameter = 1.0
        """

        pass
        #dpmm = DPMM(self.data, verbose=True)   # TODO: store the data for plots?
        #dpmm.run_sampling(max_iteration)
        #return dpmm.max_s

    def _gmm_bic(self, n_clusters):
        gmm = GaussianMixture(n_components=n_clusters, covariance_type='full')
        gmm.fit(self.data)
        return gmm.bic(self.data)

    def gmm_bic(self, max_n_clusters=100):
        """
        Gaussian Mixture Model with model selection by BIC.
        """

        gmm_bics = [self._gmm_bic(i) for i in range(1, max_n_clusters + 1)]
        logger.debug(gmm_bics)
        gmm = GaussianMixture(n_components=np.argmin(gmm_bics) + 1, covariance_type='full')
        gmm.fit(self.data)
        return gmm.predict(self.data)


class ClusteringVarMat(ClusteringNumeric):
    """
    Overview:
          self.vmatrix     ---(any clustering method)--->     self.assignment
        (variant matrix)    [ward, gmm, split_dis, dpmm]    (clustering result)

    Data format:
       self.vmatrix.shape = (self.N, self.L)
       self.N = # of units
       self.L = # of variant sites

       self.assignment.shape = (self.N)
       self.assignment[i] = Cluster index to which unit i belongs
    """

    def __init__(self, vmatrix_fname, names=None):
        data = self.load_vmatrix(vmatrix_fname)
        super().__init__(data, names)

    def load_vmatrix(self, in_fname):
        N = int(run_command(f"awk 'NR == 1 {{print length($1)}}' {in_fname}"))
        L = int(run_command(f"cat {in_fname} | wc -l"))

        # TODO: why last column of vmatrix is always 0?

        vmatrix = np.zeros((L, N), dtype=int)
        with open(in_fname, 'r') as f:
            for i, line in enumerate(f):
                vmatrix[i, :] = list(map(int, list(line.strip())))

        # Change the matrix form as (UNITS, VARIANTS)
        return vmatrix.T[1:-1, :]

    def calc_dist_mat(self):
        super().calc_dist_mat(metric="hamming")

    def generate_consensus(self, how="representative"):
        """
        Compute consensus for each cluster.

        <how> == "representative": majority vote,
              == "centroid": average (mainly for debug)
        """

        assert how in ("representative", "centroid"), "<how> must be 'representative' or 'centroid'"

        ret = {}
        index = 0
        for cluster_id, cluster_data in super().clusters():
            n_cluster = cluster_data.shape[0]
            ret[index] = (cluster_id,
                          n_cluster,
                          np.where(np.sum(cluster_data, axis=0) >= n_cluster // 2, 1, 0) if how == "representative"
                          else np.sum(cluster_data, axis=0) / n_cluster)

        return pd.DataFrame.from_dict(ret, orient="index",
                                      columns=("cluster_id", "cluster_size", "representative"))
