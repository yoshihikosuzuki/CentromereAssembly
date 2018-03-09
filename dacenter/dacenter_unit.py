import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn import cluster
from sklearn import mixture
from sklearn import decomposition

#from dpmm import *
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

    def generate_representative_units(self):
        self.r_units = []

        for i in range(1, max(self.assignment) + 1):
            units = self.vmatrix[np.where(self.assignment == i)[0], :]
            cluster_size = len(units[:, 0])
            units_sum = np.sum(units, axis=0)
            self.r_units.append(np.array([1 if x >= cluster_size / 2
                                          else 0
                                          for x in units_sum]))

    def output_representative_units(self, out_fname):
        with open(out_fname, 'w') as f:
            for r_id, r_seq in enumerate(self.r_units):
                header = "representative/%d/0_%d" % (r_id, len(r_seq))
                f.write(">%s\n%s\n" % (header, r_seq))
