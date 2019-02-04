import random
from logzero import logger
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
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go
from BITS.run import run_edlib
from BITS.utils import run_command, print_log, NoDaemonPool
import consed
from .dpmm import DPMM, DPMMCluster
#from .dpmm_oldname import Clustering, Cluster


class Clustering:
    """
    Root class of clusterings of
      1. DNA sequences (regardless of cyclic alignment and/or revcomp)
      2. variant matrix (by Consed)

    <self.assignment> is the final result of the clustering.

    Both <self.s_dist_mat> (squre distance matrix) and <self.c_dist_mat> (condensed distance matrix)
    must be set when the distance matrix is computed.
    """

    def __init__(self, input_data):
        self.data = input_data
        self.N = input_data.shape[0]   # num of data. NOTE: <input_data> must be (data, features)
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

    def plot_heatmap(self, cmap="Blues_r", figsize=(12, 12), out_fname=None):
        """
        Draw a heatmap of the (squared) distance matrix.
        """

        plt.figure(figsize=figsize)
        sns.heatmap(self.s_dist_mat, square=True, vmin=0, vmax=1, cmap=cmap)
        if out_fname is None:
            plt.show()
        else:
            plt.savefig(out_fname)

    def plot_tsne(self, coloring="sequential", figsize=(12, 12), out_fname=None):
        """
        Embed data into a two dimensional space using t-SNE.
        <self.s_dist_mat> must be precomputed.
        """

        assert hasattr(self, "s_dist_mat"), "No square distance matrix"
        assert coloring in ("sequential", "random"), "<color> must be 'sequential' or 'random'"

        trace = go.Scatter(x=self.coord[:, 0].T,
                           y=self.coord[:, 1].T,
                           text=self.assignment,
                           mode="markers",
                           marker=dict(size=3,
                                       color=self.assignment,
                                       colorscale="Rainbow",
                                       showscale=False))
        layout = go.Layout(width=700,
                           height=700,
                           hovermode='closest')
        py.iplot(go.Figure(data=[trace], layout=layout))

        """
        if not hasattr(self, "coord"):
            self.coord = TSNE(n_components=2, metric='precomputed').fit_transform(self.s_dist_mat)

        if coloring == "sequential":
            cmap = plt.get_cmap("gist_ncar")
            cols = cmap(np.array(range(self.n_clusters)) / self.n_clusters)
        plt.figure(figsize=figsize)
        for i, data in enumerate(self.clusters(return_where=True)):
            cluster_id, where = data
            plt.scatter(*self.coord[where].T, s=5,
                        c=cols[i] if coloring == "sequential" else f"#{random.randint(0, 0xFFFFFF):06x}",
                        label=f"{cluster_id} ({where.shape[0]} seqs)")
        plt.legend(fontsize=12, markerscale=4, bbox_to_anchor=(1.05, 1))

        if out_fname is None:
            plt.show()
        else:
            plt.savefig(out_fname)
        """


def __calc_dist_array(i, data):
    """
    This is just for parallelization of distance matrix calculation
    """

    # row <i> in the distance matrix
    dist_array = np.array([run_edlib(data[i],
                                     data[j],
                                     "global",
                                     cyclic=True,
                                     rc=True,
                                     only_diff=True)
                           for j in range(i + 1, data.shape[0])],
                          dtype='float32')

    logger.debug(f"Finished @ row {i}")
    return (i, dist_array)


def _calc_dist_array(args):
    s, t, data = args
    logger.debug(f"Starting row {s}-{t}")
    return [__calc_dist_array(i, data) for i in range(s, t)]


class ClusteringSeqs(Clustering):
    """
    Perform clustering of DNA sequences given.
    Both greedy approach based on cluster diameter and hierarchical clustering approach are available.

    Run <self.cluster_greedy> for former, and <self.cluster_hierarchical> for latter.
    """

    def __init__(self, input_data, cyclic=True, rc=True):
        super().__init__(input_data)
        self.cyclic = cyclic   # do cyclic alignment between two sequences
        self.rc = rc   # allow reverse complement when taking alignment

    @print_log("distance matrix calculation")
    def calc_dist_mat(self, n_core=1):
        """
        Calculate all-vs-all distance matrix between the sequences.
        Both cyclic alignment and reverse complement sequence are considered.
        Computation is performed in parallel, and the unit of parallelization
        is each row of the (triangular) distance matrix.
        """

        # Allocate a square distance matrix
        self.s_dist_mat = np.zeros((self.N, self.N), dtype='float32')

        # Split jobs while considering that weight is different for each row
        tasks = []
        n_sub = -(-((self.N - 1) * (self.N - 1)) // n_core / 2)
        s = 0
        total = 0
        for t in range(self.N - 1):
            total += self.N - 1 - t
            if total >= n_sub:
                tasks.append((s, t + 1, self.data))
                s = t + 1
                total = 0
        if total > 0:
            tasks.append((s, self.N - 1, self.data))

        with NoDaemonPool(n_core) as pool:
            for ret in pool.map(_calc_dist_array, tasks):
                for r in ret:
                    i, dist_array = r
                    self.s_dist_mat[i, i + 1:] = self.s_dist_mat[i + 1:, i] = dist_array

        # Generate a condensed matrix
        self.c_dist_mat = squareform(self.s_dist_mat)

    def generate_consensus(self):
        """
        Determine a representative sequence for each cluster.
        """

        ret = {}
        index = 0
        for cluster_id, seqs in super().clusters():
            cons_seq = consed.consensus([seq if i == 0
                                         else run_edlib(seqs.iloc[0],
                                                        seq,
                                                        "global",
                                                        cyclic=True,
                                                        rc=True,
                                                        return_seq=True).seq
                                         for i, seq in enumerate(seqs)],
                                        n_iter=3)
            ret[index] = (cluster_id, seqs.shape[0], len(cons_seq), cons_seq)
            index += 1

        return pd.DataFrame.from_dict(ret, orient="index",
                                      columns=("cluster_id", "cluster_size", "length", "sequence"))


class ClusteringNumeric(Clustering):
    """
    A child class which additionally has clustering methods only for numerical data.
    """

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

        dpmm = DPMM(self.data, verbose=True)   # TODO: store the data for plots?
        dpmm.run_sampling(max_iteration)
        return dpmm.max_s

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

    def __init__(self, vmatrix_fname):
        super().__init__(self.load_vmatrix(vmatrix_fname))

    def load_vmatrix(self, in_fname):
        N = int(run_command(f"awk 'NR == 1 {{print length($1)}}' {in_fname}"))
        L = int(run_command(f"cat {in_fname} | wc -l"))

        # TODO: why last column of vmatrix is always 0?

        vmatrix = np.zeros((L, N), dtype=int)
        with open(in_fname, 'r') as f:
            for i, line in enumerate(f):
                vmatrix[i, :] = list(map(int, list(line.strip())))

        # Change the matrix form as (UNITS, VARIANTS)
        return vmatrix.T

    def calc_dist_mat(self):
        # Used for hierarchical clustering and/or t-SNE plot
        self.c_dist_mat = pdist(self.data, metric='hamming')
        self.s_dist_mat = squareform(self.c_dist_mat)

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
