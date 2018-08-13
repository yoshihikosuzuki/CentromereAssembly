import argparse
import os
import pickle
from logzero import logger
from collections import Counter
from multiprocessing import Pool
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.cluster import KMeans, Birch
from sklearn.mixture import GaussianMixture
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go

from BITS.utils import run_command, revcomp
from BITS.core import run_edlib

from .dpmm import DPMM, DPMMCluster
#from .dpmm_oldname import Clustering, Cluster


class Clustering:
    """
    Super class of clusterings of
      1. DNA sequences (regardless of cyclic alignment and/or revcomp)
      2. distance matrix (distance = sequence difference)
      3. variant matrix (by Consed)
    """

    def __init__(self, input_data):
        self.data = input_data
        self.N = input_data.shape[0]   # NOTE: this must be data size!
        self.assignment = np.full(self.N, -1, dtype='int8')   # cluster assignment for each data

    def calc_clusters(self):
        """
        From <self.assignment>, generate the clusters with data belonging to one of them.
        """

        self.clusters = {}
        for cluster_id in set(self.assignment):
            self.clusters[cluster_id] = self.data[np.where(self.assignment == cluster_id)]


class ClusteringSeqs(Clustering):
    def __init__(self, input_data, cyclic=True, revcomp=True):
        super().__init__(input_data)
        self.cyclic = cyclic   # do cyclic alignment between two sequences
        self.revcomp = revcomp   # allow reverse complement when taking alignment

    #def cluster_greedy(self):

    def calc_dist_array(self, i):
        logger.debug(f"Started job: row {i}, columns {i + 1}-{self.N - 1}")
        dist_array = np.empty(self.N - i - 1, dtype='float32')   # of row <i> in the dist matrix

        query = self.data[i]
        query_rc = revcomp(query)
        for index, j in enumerate(range(i + 1, self.N - 1)):
            target = self.data[j] * 2
            alignment_f = run_edlib(query, target, mode="glocal")   # TODO: no need of <task="path"> in edlib
            if alignment_f["diff"] <= 0.3:
                alignment = alignment_f   # TODO: needed here is only <diff> (maybe should check length?)
            else:
                alignment_rc = run_edlib(query_rc, target, mode="glocal")
                alignment = alignment_f if alignment_f["diff"] <= alignment_rc["diff"] else alignment_rc
            dist_array[index] = alignment["diff"]

        logger.debug(f"Ended job: row {i}")

    def calc_dist_mat(self, n_core=1):
        """
        Calculate all-vs-all distance matrix between the sequences.
        Cyclic alignment, allowing reverse complement of one sequence during alignment,
        and parallelization are supported.
        """

        logger.info(f"Starting calculation of distance matrix with {n_core} cores")

        exe_pool = Pool(n_core)
        self.dist_matrix = np.zeros((self.N, self.N), dtype='float32')

        # NOTE: each row in the matrix is unit of the parallelized tasks
        logger.debug("Scattering tasks")

        for ret in exe_pool.imap(self.calc_dist_array, np.arange(self.N - 1)):
            i, dist_array = ret
            self.dist_matrix[i, i + 1:] = self.dist_matrix[i + 1:, i] = dist_array
            logger.debug(f"Finished row {i}")

        exe_pool.close()
        logger.debug("Finished all tasks")

    #def cluster_hierarchical(self, method="ward"):


#class ClusteringDistMat:

class ClusteringVarMat:
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
        gmm = mixture.GaussianMixture(n_components=n_components,
                                      covariance_type='full')
        gmm.fit(self.vmatrix)
        return gmm.bic(self.vmatrix)

    # Gaussian Mixture Model with model selection by BIC
    def gmm(self, max_components=100):
        gmm_bics = [self.gmm_bic(i) for i in range(1, max_components + 1)]
        print(gmm_bics)
        gmm = mixture.GaussianMixture(n_components=np.argmin(gmm_bics) + 1,
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
        self.assignment = cluster.KMeans(n_clusters=n_components).fit_predict(self.vmatrix)

    def gmm_split(self, n_components):
        gmm = mixture.GaussianMixture(n_components=n_components,
                                      covariance_type='full')
        gmm.fit(self.vmatrix)
        self.assignment = gmm.predict(self.vmatrix)

    def birch(self, n_components):
        self.assignment = cluster.Birch(n_clusters=n_components).fit_predict(self.vmatrix)

    def nmf(self, n_components):
        nmf = decomposition.NMF(n_components=n_components)
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

    def save_pickle(self, out_fname):
        with open(out_fname, 'wb') as f:
            pickle.dump(self, f)


def run_consed(in_fname):
    out_fname_prefix = os.path.splitext(in_fname)[0]

    print("[INFO] If Consed failed with \"Cannot align to first sequence\", "
          "try \"$ sed -i -e \"<read_id + 1>d\" <aligned_units_fname>\".")

    # Consed output with variant matrix
    command = "consed -V -w1000000 %s" % in_fname
    consed_fname = out_fname_prefix + ".consed"
    with open(consed_fname, 'w')as f:
        f.write(run_command(command))

    # Extract consensus sequence
    command = ("awk 'BEGIN {seq = \"\"} $0 == \"\" {exit} {seq = seq $0} END {print seq}' %s"
               % consed_fname)
    cons_seq = run_command(command).strip()
    with open(out_fname_prefix + ".consensus.fasta", 'w') as f:
        f.write(">consensus/0/0_%d\n%s\n" % (len(cons_seq), cons_seq))

    # Extract variant matrix
    command = "awk 'NR > 5 {print $NF}' %s" % consed_fname
    with open(out_fname_prefix + ".V", 'w') as f:
        f.write(run_command(command))

    return out_fname_prefix + ".V"


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
