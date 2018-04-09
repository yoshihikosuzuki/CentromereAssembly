import pickle
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA, KernelPCA
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go
from collections import Counter

from dacenter_run_peak import Peaks
from dacenter_run_unit import Clustering
from dpmm import DPMM, DPMMCluster
#from dpmm_oldname import Clustering, Cluster

plt.style.use('ggplot')
py.init_notebook_mode()


class PeaksViewer:
    def __init__(self, unit_fasta_fname):
        self.peaks = Peaks(None, None)
        self.peaks.load_units(unit_fasta_fname)

    def hist_units(self, min_len=50, max_len=450, band_width=10):
        unit_lens = np.array([len(x[1]) for x in self.peaks.units])

        data = [go.Histogram(x=unit_lens,
                             xbins=dict(start=min_len, end=max_len, size=1))]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Frequency"))
        py.iplot(go.Figure(data=data, layout=layout))

        ul = unit_lens[(min_len < unit_lens) & (unit_lens < max_len)]
        ls = np.linspace(min_len, max_len, max_len - min_len + 1, dtype=int)

        kde = KernelDensity(kernel='gaussian', bandwidth=band_width).fit(ul.reshape(-1, 1))
        dens = np.exp(kde.score_samples(ls.reshape(-1, 1)))

        data = [go.Scatter(x=ls, y=dens)]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Density"))
        py.iplot(go.Figure(data=data, layout=layout))


class ClusteringViewer:
    def __init__(self, in_pkl_fname):
        with open(in_pkl_fname, 'rb') as f:
            self.clustering = pickle.load(f)

        # NOTE: Viewer needs clustering has at least:
        #    assignment (max_s in dpmm.oldname.py)
        #    vmatrix (x in dpmm.oldname.py)
        #    N, L
        #    r_units

    # Transition of posterior probability and # of clusters
    # Can be used only when DPMM clustering was performed
    def plot_transition(self, start=0, end=None, figsize=(10, 8)):
        if not hasattr(self.clustering, "dpmm"):
            return

        fig, ax1 = plt.subplots(figsize=figsize)
        ax1.plot([x[0] for x in self.clustering.dpmm.posterior_ncluster[start:end]], color='r')
        ax2 = ax1.twinx()
        ax2.plot([x[1] for x in self.clustering.dpmm.posterior_ncluster[start:end]], color='b')
        plt.grid(False)
        plt.show()

    def hist_cluster_size(self, bins=100, n_print=5):
        cluster_size = Counter(self.clustering.assignment)
        print("%d largest clusters:" % n_print)
        for i, pair in enumerate(sorted(cluster_size.items(),
                                        key=lambda x: x[1],
                                        reverse=True)):
            if i >= n_print:
                break
            print("cluster %d (%d units)" % (pair[0], pair[1]))
        print("%d smallest clusters:" % n_print)
        for i, pair in enumerate(sorted(cluster_size.items(),
                                        key=lambda x: x[1])):
            if i >= n_print:
                break
            print("cluster %d (%d units)" % (pair[0], pair[1]))
        plt.hist(list(cluster_size.values()), bins=bins)
        plt.show()

    # Single heatmap with partitions of the clusters
    # Only for relatively small (<1000) dataset
    def heatmap_partition(self, partition_width=3, figsize=(18, 15)):
        if self.clustering.N > 2000:
            print("[ERROR] Data size is too large to draw single heatmap. "
                  "Use heatmap_clusters() instead.")
            return

        d = np.full((self.clustering.L,
                     self.clustering.N + partition_width * (len(set(self.clustering.assignment)) + 1)),
                    -1,
                    dtype=int)
        start = partition_width
        for idx in set(self.clustering.assignment):
            cluster_data = self.clustering.vmatrix.T[:, np.where(self.clustering.assignment == idx)[0]]
            d[:, start:start + cluster_data.shape[1]] = cluster_data
            start += cluster_data.shape[1] + partition_width

        plt.subplots(figsize=figsize)
        sns.heatmap(d, cbar=False)
        plt.show()

    # Heatmap of a single cluster specified by c_idx
    def heatmap_cluster(self, c_idx, figsize=(8, 6)):
        plt.subplots(figsize=figsize)
        sns.heatmap(self.clustering.vmatrix[np.where(self.clustering.assignment == c_idx)[0], :], cbar=False)
        plt.show()

    # Distinct heatmaps of all clusters
    def heatmap_clusters(self, figsize=(8, 6)):
        for c_idx in set(self.clustering.assignment):
            self.heatmap_cluster(c_idx, figsize)

    # Dendrogram of the representative units by Hamming distance and Ward's method
    def dendrogram_representatives(self):
        hc_result = linkage(pdist(list(self.clustering.r_units.values()), metric='hamming'), method='ward')
        plt.figure(figsize=(18, 10))
        dendrogram(hc_result)
        plt.show()

    # Scatter plot of raw units as 2D Hamming kernel PCA
    def scatter_units(self, figsize=(12, 12), only_representatives=False):
        # Add representative units into the variant matrix
        vm = self.clustering.vmatrix.copy()
        rm = np.zeros((len(self.clustering.r_units), self.clustering.L), dtype=int)
        for i, r_unit in enumerate(sorted(self.clustering.r_units.items(), key=lambda x: x[0])):
            rm[i] = r_unit[1]   # add in ascending order of cluster index
        data = np.concatenate((vm, rm), axis=0)

        coord = KernelPCA(n_components=2, kernel='precomputed').fit_transform(1 - squareform(pdist(data, metric='hamming')))

        # NOTE: normal PCA gives almost same result because of sparsity
        #coord = PCA(n_components=2).fit_transform(data)

        fig, ax = plt.subplots(figsize=figsize)
        for i, c_idx in enumerate(sorted(list(set(self.clustering.assignment)))):
            col = "#%06X" % (int(0xFFFFFF / (len(set(self.clustering.assignment)) - 1) * i))
            if not only_representatives:
                ax.scatter(coord[np.where(self.clustering.assignment == c_idx)[0], 0], coord[np.where(self.clustering.assignment == c_idx)[0], 1], s=20, c=col, label=str(c_idx), marker=".")
            marker = "$%s$" % c_idx if only_representatives else "*"
            ax.scatter(coord[-(len(self.clustering.r_units) - i), 0], coord[-(len(self.clustering.r_units) - i), 1], s=500, c=col, marker=marker)   # representative units
        if not only_representatives:
            ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1.0), prop={"size": 12})
        plt.show()
