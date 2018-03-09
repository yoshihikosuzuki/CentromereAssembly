import argparse
import os
import pickle
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn import cluster
from sklearn import mixture
from sklearn import decomposition

from dpmm import DPMM, DPMMCluster
from datruf_utils import run_command


class Clustering:
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

    # Hierarchical, Ward, 0.7
    def ward(self, threshold=0.7):
        hc_ward = linkage(pdist(self.vmatrix, metric='hamming'), method='ward')
        self.assignment = np.array(fcluster(hc_ward, threshold, criterion='distance'))

    # DPMM, alpha=1.0
    def dpmm(self, max_iteration=1000000):   # TODO: tune parameters (max_iteration, converge_threshold_*)
        self.dpmm = DPMM(self.vmatrix)
        self.dpmm.run_sampling(max_iteration)
        self.assignment = self.dpmm.max_s

    def generate_representative_units(self):
        self.r_units = []

        for c_idx in set(self.assignment):
            cluster_units = self.vmatrix[np.where(self.assignment == c_idx)[0], :]
            cluster_size = cluster_units.shape[0]
            units_sum = np.sum(cluster_units, axis=0)
            self.r_units.append(np.array([1 if x >= cluster_size / 2
                                          else 0
                                          for x in units_sum]))

    def output_representative_units(self, out_fname):
        with open(out_fname, 'w') as f:
            for r_id, r_seq in enumerate(self.r_units):
                header = "representative/%d/0_%d" % (r_id, len(r_seq))
                f.write(">%s\n%s\n" % (header, r_seq))

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
    #clustering.ward()
    clustering.dpmm()
    print("[INFO] %d clusters were generated." % len(set(clustering.assignment)))
    clustering.generate_representative_units()
    #clustering.output_representative_units("repr.ward.fasta")
    #clustering.save_pickle("ward.pkl")
    clustering.output_representative_units("repr.dpmm.fasta")
    clustering.save_pickle("dpmm.pkl")


def load_args():
    parser = argparse.ArgumentParser(
        description=("Run dacenter for clustering of units. You must specify "
                     "one of the following arguments: [--aligned_units_fname] "
                     "or [--variant_matrix_fname]."))

    parser.add_argument(
        "--aligned_units_fname",
        default=None,
        help=(".seq file of the (start-aligned) unit sequences [None]"))

    parser.add_argument(
        "--variant_matrix_fname",
        default=None,
        help=("Precomputed variant matrix file [None]"))

    args = parser.parse_args()

    if args.aligned_units_fname is None and args.variant_matrix_fname is None:
        print("[ERROR] You must specify one of the following arguments:\n"
              "[--aligned_units_fname] or [--variant_matrix_fname].")
        exit(1)

    return args


if __name__ == "__main__":
    main()
