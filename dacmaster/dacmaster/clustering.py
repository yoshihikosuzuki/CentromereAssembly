import sys
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
import plotly.offline as py
import plotly.graph_objs as go
from BITS.seq import revcomp
from BITS.run import run_edlib
from BITS.utils import run_command, print_log, NoDaemonPool
import consed
from .dpmm import DPMM, DPMMCluster
#from .dpmm_oldname import Clustering, Cluster


class Clustering:
    """
    Super class of clusterings of
      1. DNA sequences (regardless of cyclic alignment and/or revcomp)
      2. variant matrix (by Consed)

    <self.assignment> is the final result of the clustering.
    """

    def __init__(self, input_data):
        self.data = input_data
        self.N = input_data.shape[0]   # NOTE: <input_data> must be this kind of data
        self.assignment = np.full(self.N, -1, dtype='int8')   # cluster assignment for each data
        self.hc_result_precomputed = {}   # used to avoid re-calculation of hierarchical clustering

    def cluster_hierarchical(self, method, criterion, threshold):
        """
        <method> = {single, complete, average, weighted, centroid, median, ward (default)}
        <criterion> = {inconsistent, distance (default), maxclust, monocrit, maxclust_monocrit}
        <threshold> = threshold in distance criterion

        Before calling this, <self.hc_input> must be set with a condensed distance matrix.

        NOTE: Currently this supports only distance criterion.
        """

        # Since computation of linkage is heavy to some extent, save it for the next time
        if method in self.hc_result_precomputed:
            self.hc_result = self.hc_result_precomputed[method]
        else:
            self.hc_result = linkage(self.hc_input, method=method)
            self.hc_result_precomputed[method] = self.hc_result

        # Calculate cluster assignment for each data
        self.assignment = np.array(fcluster(self.hc_result, t=threshold, criterion=criterion))

    def dendrogram(self, method):
        """
        Show dendrogram of data with <method>.
        """

        if method in self.hc_result_precomputed:
            self.hc_result = self.hc_result_precomputed[method]
        else:
            self.hc_result = linkage(self.hc_input, method=method)
            self.hc_result_precomputed[method] = self.hc_result

        plt.figure(figsize=(18, 10))
        dendrogram(self.hc_result)
        plt.show()

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
        Be careful it is a generater.
        """

        for cluster_id in sorted(list(set(self.assignment))):
            where = np.where(self.assignment == cluster_id)[0]
            if return_where:
                yield (cluster_id, where)
            else:
                yield (cluster_id, self.data[where])

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
    def calc_dist_mat(self, n_core):
        """
        Calculate all-vs-all distance matrix between the sequences.
        Both cyclic alignment and reverse complement sequence are considered.
        Computation is performed in parallel, and the unit of parallelization
        is each row of the (triangular) distance matrix.
        """

        self.dist_matrix = np.zeros((self.N, self.N), dtype='float32')

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

        exe_pool = NoDaemonPool(n_core)
        for ret in exe_pool.map(_calc_dist_array, tasks):
            for r in ret:
                i, dist_array = r
                self.dist_matrix[i, i + 1:] = self.dist_matrix[i + 1:, i] = dist_array
        exe_pool.close()

    def cluster_hierarchical(self, method="ward", criterion="distance", threshold=0.7, n_core=1):
        if not hasattr(self, "dist_matrix"):
            self.calc_dist_mat(n_core)

        if not hasattr(self, "hc_input"):
            self.hc_input = squareform(self.dist_matrix)

        super().cluster_hierarchical(method, criterion, threshold)

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
                                                        return_seq=True)["seq"]
                                         for i, seq in enumerate(seqs)],
                                        n_iter=3)
            ret[index] = (cluster_id, seqs.shape[0], len(cons_seq), cons_seq)
            index += 1

        return pd.DataFrame.from_dict(ret, orient="index",
                                      columns=("cluster_id", "cluster_size", "length", "sequence"))

    def dendrogram(self, method="ward"):
        super().dendrogram(method)

    def plot_tsne(self, figsize=(12, 12), leg_marker_size=200, out_fname=None):
        if not hasattr(self, "coord"):
            self.coord = TSNE(n_components=2, metric='precomputed').fit_transform(self.dist_matrix)

        fig, ax = plt.subplots(figsize=figsize)
        for i, data in enumerate(super().clusters(return_where=True)):
            cluster_id, where = data
            ax.scatter(self.coord[where, 0], self.coord[where, 1], marker=".", s=20,
                       c=f"#{random.randint(0, 0xFFFFFF):06x}",
                       label=f"{cluster_id} ({where.shape[0]} seqs)")
        ax.legend(loc="upper right", bbox_to_anchor=(1.2, 1.0), prop={"size": 12})
        leg = ax.get_legend()
        for handle in leg.legendHandles:
            handle._sizes = [leg_marker_size]

        if out_fname is not None:
            plt.savefig(out_fname)
        else:
            plt.show()

    # ------------------------------
    # Greedy clustering from here...
    # ------------------------------

    def propose_cluster(self, remaining_read_id, cluster_id, out_dir):
        """
        From the remaining consensus units, first randomly choose a unit and
        then collect all units within diameter of 10% sequence difference from
        it.
        """

        cluster_units = np.zeros(self.units_consensus.shape[0], dtype=int)   # <cluster_id> + 1 for cluster members, otherwise 0
        read_id = random.choice(remaining_read_id)
        query = self.units_consensus["sequence"].iloc[read_id]
        
        out_fname = out_dir +  "/input.seq"   # temporary file for Consed input
        seq_writer = open(out_fname, 'w')
    
        # Output the first sequence as seed for consensus
        seq_writer.write(f"{query}\n")
    
        for j in remaining_read_id:
            target = self.units_consensus["sequence"].iloc[j]
            
            # forward
            seq_f = target * 2
            alignment_f = run_edlib(query, seq_f, mode="glocal")
    
            # reverse complement
            seq_rc = revcomp(seq_f)
            alignment_rc = run_edlib(query, seq_rc, mode="glocal")
    
            if alignment_f["diff"] <= alignment_rc["diff"]:
                alignment = alignment_f
                seq = seq_f
            else:
                alignment = alignment_rc
                seq = seq_rc
    
            if alignment["diff"] < 0.1:   # 0.05とかにしてもいいかも
                cluster_units[j] = cluster_id + 1
                seq = seq[alignment["start"]:alignment["end"]]   # mapped area
                seq_writer.write(f"{seq}\n")
    
        seq_writer.close()
        return (cluster_units,
                run_command(f"consed {out_fname}").replace('\n', ''))

    def fix_cluster_assignment(self, remaining_read_id, cluster_id, consensus_seq):   # TODO: can I merge this and above one into a single method?
        """
        Given a set of consensus units proposed, return an actual cluster
        which has diameter of 10% sequence difference from consensus of the
        consensus units. The cluster size might become very smaller than the
        proposed one.
        """

        cluster_units = np.zeros(self.units_consensus.shape[0], dtype=int)   # 1 for cluster members, otherwise 0
        read_id = random.choice(remaining_read_id)
        query = self.units_consensus["sequence"].iloc[read_id]
    
        for j in remaining_read_id:
            target = self.units_consensus["sequence"].iloc[j]
    
            # forward
            seq_f = target * 2
            alignment_f = run_edlib(query, seq_f, mode="glocal")
    
            if alignment_f["diff"] < 0.1:
                cluster_units[j] = cluster_id + 1   # add 1 to neutralize the default value -1 in the <assignment>
            else:
                # reverse complement
                seq_rc = revcomp(seq_f)
                alignment_rc = run_edlib(query, seq_rc, mode="glocal")
                if alignment_rc["diff"] < 0.1:
                    cluster_units[j] = cluster_id + 1
                    
        return cluster_units

    def cluster_greedy(self):   # TODO: refactor so that this works
        # Clustering of the consensus untis
        # Since what we are doing for now is just determine rough representative monomers,
        # this task is rather just to eliminate redundant consensus units.
        self.assignment = np.full(self.units_consensus.shape[0], -1, dtype=int)   # cluster id assignment for each consensus unit
        self.clusters = {}   # mater unit sequences
        cluster_id = 0
        remaining_read_id = np.arange(self.units_consensus.shape[0])   # array of indices of consensus units still remaining

        while True:
            logger.info(f"propose {cluster_id}")
            cluster_units, consensus_seq = self.propose_cluster(remaining_read_id,
                                                                cluster_id,
                                                                f"tmp/")
            
            # XXX: consensus sequence may be "EXITING", "** WARNING", or Nothing!!!!
            # TODO: maybe I should first debug Consed
            
            cluster_size = cluster_units[cluster_units > 0].shape[0]
            print(cluster_id, cluster_size)
            
            #print(consensus_seq)
            if cluster_size >= self.units_consensus.shape[0] * 0.01:
                print("accepted")
                
                # Re-align the consensus sequence to the remaining consensus units
                cluster_units = self.fix_cluster_assignment(remaining_read_id,
                                                            cluster_id,
                                                            consensus_seq)   # TODO: more efficient way?
                
                cluster_size = cluster_units[cluster_units > 0].shape[0]
                print("actual cluster size =", cluster_size)
                
                if cluster_size >= self.units_consensus.shape[0] * 0.01:   # filter by cluster size again
                    print("accepted")
                    self.clusters[cluster_id] = consensus_seq
                    self.assignment += cluster_units
                    cluster_id += 1
                else:
                    print("rejected")
            else:
                print("rejected")
                self.assignment -= cluster_units
                
            remaining_read_id = np.where(self.assignment == -1)[0]
            print("# of remaining units =", remaining_read_id.shape[0])
            #print(clusters)
            if remaining_read_id.shape[0] < self.units_consensus.shape[0] * 0.1:
                break
        
        self.assignment[np.where(self.assignment < 0)] = -1

        logger.info(f"{len(self.clusters)} representative unit candidates:\n{self.clusters}")
        # remove irregular sequences
        """
        del_ids = set()
        for cluster_id, cluster_unit in self.clusters.items():
            if len(cluster_unit) == 0 or cluster_unit[0] == 'W' or cluster_unit[0] == '*':
                del_ids.add(cluster_id)
        for cluster_id in del_ids:
            del self.clusters[cluster_id]
        """
        
        self.clusters = pd.DataFrame.from_dict(self.clusters, orient="index", columns=["sequence"])
        for i in range(self.clusters.shape[0]):
            if len(self.clusters.loc[i, "sequence"]) == 0 or self.clusters.loc[i, "sequence"][0] == '*' or self.clusters.loc[i, "sequence"] == 'W':
                self.clusters = self.clusters.drop(i)
        self.clusters = self.clusters.reset_index(drop=True)

        if self.clusters.shape[0] == 0:
            logger.info("No representatives found. Exit.")
            sys.exit(1)

        print(self.clusters)
        # TODO: add cluster size column and sort by it

        # Remove redundant untis which are similar to another one
        # After that, the remaining centroids are representative units
        dm = np.zeros((self.clusters.shape[0], self.clusters.shape[0]), dtype=float)
        for i in range(self.clusters.shape[0] - 1):
            for j in range(1, self.clusters.shape[0]):
                #print(f"# {i} vs {j} (f)")
                alignment = run_edlib(self.clusters.loc[i, "sequence"], self.clusters.loc[j, "sequence"] * 2, mode='glocal')
                f = alignment["diff"]
                #print(f"# {i} vs {j} (rc)")
                alignment = run_edlib(self.clusters.loc[i, "sequence"], revcomp(self.clusters.loc[j, "sequence"] * 2), mode='glocal')
                rc = alignment["diff"]
                
                dm[i, j] = dm[j, i] = min([f, rc])

        self.representative_units = {self.clusters.loc[0, "sequence"]}
        for i in range(1, self.clusters.shape[0]):
            flag_add = 1
            for j in range(i):
                if dm[i, j] < 0.1:
                    flag_add = 0
                    break
            if flag_add == 1:
                self.representative_units.add(self.clusters.loc[i, "sequence"])

        logger.info(f"{len(self.representative_units)} after redundancy removal:\n{self.representative_units}")


class ClusteringVarMat:   # TOOD: change to a child class of Clustering
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
        self.load_vmatrix(vmatrix_fname)

    def load_vmatrix(self, in_fname):
        self.N = int(run_command("awk 'NR == 1 {print length($1)}' %s"
                                 % in_fname))
        self.L = int(run_command("cat %s | wc -l" % in_fname))

        # TODO: last column of vmatrix is always 0 (i.e. consensus seuqence?). should remove it?

        self.vmatrix = np.zeros((self.L, self.N), dtype=int)
        with open(in_fname, 'r') as f:
            for i, line in enumerate(f):
                self.vmatrix[i, :] = list(map(int, list(line.strip())))

        # Keep the matrix as the [UNITS, VARIANTS] form
        self.vmatrix = self.vmatrix.T

    # Run clustering not requiring specification of # of clusters
    def run_clustering(self, method):
        if method == "dpmm":
            self.dpmm()
        elif method == "ward":
            self.ward()
        elif method == "gmm":
            self.gmm()

    # Dirichlet Process (Multinomial) Mixture Model
    # concentration hyperparameter = 1.0
    def dpmm(self, max_iteration=1000000):   # TODO: tune parameters (max_iteration, converge_threshold_*)
        self.dpmm = DPMM(self.vmatrix, verbose=True)
        self.dpmm.run_sampling(max_iteration)
        self.assignment = self.dpmm.max_s

    # Hierarchical clustering, Ward method
    # distance threshold = 0.7
    def ward(self, distance_threshold=0.7):
        hc_ward = linkage(pdist(self.vmatrix,
                                metric='hamming'),
                          method='ward')
        self.assignment = np.array(fcluster(hc_ward,
                                            distance_threshold,
                                            criterion='distance'))

    def gmm_bic(self, n_components):
        gmm = GaussianMixture(n_components=n_components,
                              covariance_type='full')
        gmm.fit(self.vmatrix)
        return gmm.bic(self.vmatrix)

    # Gaussian Mixture Model with model selection by BIC
    def gmm(self, max_components=100):
        gmm_bics = [self.gmm_bic(i) for i in range(1, max_components + 1)]
        print(gmm_bics)
        gmm = GaussianMixture(n_components=np.argmin(gmm_bics) + 1,
                              covariance_type='full')
        gmm.fit(self.vmatrix)
        self.assignment = gmm.predict(self.vmatrix)

    # Run split (i.e. 1 cluster -> 2 clusters)
    def run_split(self, method, n_components=2):
        if method == "kmeans":
            self.kmeans(n_components)
        elif method == "gmm":
            self.gmm_split(n_components)
        elif method == "birch":
            self.birch(n_components)
        elif method == "nmf":
            self.nmf(n_components)

    def kmeans(self, n_components):
        self.assignment = KMeans(n_clusters=n_components).fit_predict(self.vmatrix)

    def gmm_split(self, n_components):
        gmm = GaussianMixture(n_components=n_components,
                              covariance_type='full')
        gmm.fit(self.vmatrix)
        self.assignment = gmm.predict(self.vmatrix)

    def birch(self, n_components):
        self.assignment = Birch(n_clusters=n_components).fit_predict(self.vmatrix)

    def nmf(self, n_components):
        nmf = NMF(n_components=n_components)
        W = nmf.fit_transform(self.vmatrix)
        H = nmf.components_
        print(W)
        print(H)

    def generate_representative_units(self):
        self.r_units = {}
        for c_idx in set(self.assignment):
            cluster_units = self.vmatrix[np.where(self.assignment == c_idx)[0], :]
            cluster_size = cluster_units.shape[0]
            self.r_units[c_idx] = np.where(np.sum(cluster_units, axis=0) >= cluster_size // 2, 1, 0)

    def generate_centroid_units(self):
        self.c_units = {}
        for c_idx in set(self.assignment):
            cluster_units = self.vmatrix[np.where(self.assignment == c_idx)[0], :]
            cluster_size = cluster_units.shape[0]
            self.c_units[c_idx] = np.sum(cluster_units, axis=0) / cluster_size

    def output_representative_units(self, out_fname):
        with open(out_fname, 'w') as f:
            for c_idx, r_seq in self.r_units.items():
                f.write(">repr/%d/0_%d\n%s\n" % (c_idx, len(r_seq), r_seq))


"""
def main():
    args = load_args()

    # TODO: generate variant graph as well?
    if args.aligned_units_fname is not None:
        args.variant_matrix_fname = run_consed(args.aligned_units_fname)

    clustering = Clustering(args.variant_matrix_fname)

    clustering.run_clustering(args.clustering_method)
    print("[INFO] %d clusters were generated by %s."
          % (len(set(clustering.assignment)),
             args.clustering_method))
    out_prefix = args.clustering_method   # TODO: change to appropriate one

    #clustering.run_split(args.split_method, n_components=12)
    #print("[INFO] %d clusters were generated by %s."
    #      % (len(set(clustering.assignment)),
    #         args.split_method))
    #out_prefix = args.split_method   # TODO: change to appropriate one

    clustering.generate_representative_units()
    #clustering.output_representative_units("%s.repr.fasta" % out_prefix)
    clustering.save_pickle("%s.pkl" % out_prefix)


def load_args():
    parser = argparse.ArgumentParser(
        description=("Run dacenter for clustering of units. You must specify "
                     "one of the following arguments: [-a] or [-v]."))

    parser.add_argument(
        "-u",
        "--aligned_units_fname",
        default=None,
        help=(".seq file of the (start-aligned) unit sequences [None]"))

    parser.add_argument(
        "-v",
        "--variant_matrix_fname",
        default=None,
        help=("Precomputed variant matrix file [None]"))

    parser.add_argument(
        "-c",
        "--clustering_method",
        choices=["dpmm", "ward", "gmm"],
        default="dpmm",
        help=("Clusetering method to be used [dpmm]"))

    parser.add_argument(
        "-s",
        "--split_method",
        choices=["kmeans", "gmm", "birch", "nmf"],
        default="dpmm",
        help=("Clusetering method to be used [dpmm]"))

    args = parser.parse_args()

    if args.aligned_units_fname is None and args.variant_matrix_fname is None:
        print("[ERROR] You must specify one of the following arguments:\n"
              "[-a] or [-v].")
        exit(1)

    return args


if __name__ == "__main__":
    main()
"""
