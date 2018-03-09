import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go
from collections import Counter

from datruf_utils import run_command
from dpmm import *

plt.style.use('ggplot')
py.init_notebook_mode()



    def estimate_peaks(self, min_len=50, max_len=450, band_width=10):
        # TODO: how to determine min/max_len? peak + density threshold?
    
        unit_lens = np.array([len(x[1]) for x in units])
    
        if visualize:
            data = [go.Histogram(x=unit_lens,
                                 xbins=dict(start=20, size=1))]
            layout = go.Layout(xaxis=dict(title="Unit length"),
                               yaxis=dict(title="Frequency"))
            py.iplot(go.Figure(data=data, layout=layout))
    
        ul = unit_lens[(min_len < unit_lens) & (unit_lens < max_len)]
        ls = np.linspace(min_len, max_len, max_len - min_len + 1, dtype=int)
    
        kde = KernelDensity(kernel='gaussian', bandwidth=band_width).fit(ul.reshape(-1, 1))
        dens = np.exp(kde.score_samples(ls.reshape(-1, 1)))
    
        if visualize:
            data = [go.Scatter(x=ls, y=dens)]
            layout = go.Layout(xaxis=dict(title="Unit length"),
                               yaxis=dict(title="Density"))
            py.iplot(go.Figure(data=data, layout=layout))
    
        # Naive peak detection
        peaks = []
        for i in range(1, len(dens) - 1):
            if (dens[i - 1] < dens[i]) and (dens[i] > dens[i + 1]):
                peaks.append(ls[i])
                print("[INFO] peak detected:", ls[i], dens[i])
    
        return peaks



class Viewer:
    def __init__(self, clustering_pickle):
        self.clustering = load_clustering(clustering_pickle)   # TODO: cope with other types of clustering

        # NOTE: Viewer needs clustering has at least:
        #    max_s
        #    x
        #    N, L

    def hist_cluster_size(self, bins=100, n_print=5):
        #plt.hist(self.clustering.max_s, bins=(max(self.clustering.max_s) - min(self.clustering.max_s) + 1))
        cluster_size = Counter(self.clustering.max_s)
        print("%d largest clusters:" % n_print)
        for i, pair in enumerate(sorted(cluster_size.items(), key=lambda x: x[1], reverse=True)):
            if i >= n_print:
                break
            print("cluster %d (%d units)" % (pair[0], pair[1]))
        print("%d smallest clusters:" % n_print)
        for i, pair in enumerate(sorted(cluster_size.items(), key=lambda x: x[1])):
            if i >= n_print:
                break
            print("cluster %d (%d units)" % (pair[0], pair[1]))
        plt.hist(list(cluster_size.values()), bins=bins)
        plt.show()

    # Single heatmap with partitions of the clusters
    # Only for relatively small (<1000) dataset
    def heatmap_partition(self, partition_width=3, figsize=(18, 15)):
        if self.clustering.N > 1000:   # TODO: find appropriate value
            print("[ERROR] Data size is too large to draw single heatmap. Use heatmap_clusters() instead.")
            return

        d = np.full((self.clustering.L, self.clustering.N + partition_width * (len(set(self.clustering.max_s)) + 1)), -1, dtype=int)
        start = partition_width
        for idx in set(self.clustering.max_s):
            cluster_data = self.clustering.x.T[:, np.where(self.clustering.max_s == idx)[0]]
            d[:, start:start + cluster_data.shape[1]] = cluster_data
            start += cluster_data.shape[1] + partition_width

        plt.subplots(figsize=figsize)
        sns.heatmap(d, cbar=False)
        plt.show()

    # Heatmap of a single cluster specified by c_idx
    def heatmap_cluster(self, c_idx, figsize=(8, 6)):
        plt.subplots(figsize=figsize)
        sns.heatmap(self.clustering.x[np.where(self.clustering.max_s == c_idx)[0], :], cbar=False)
        plt.show()

    # Distinct heatmaps of all clusters
    def heatmap_clusters(self, figsize=(8, 6)):
        for i in set(self.clustering.max_s):
            self.heatmap_cluster(i, figsize)

    # Dendrogram of the representative units by Hamming distance and Ward's method
    def dendrogram_representatives(self):   # TODO: separate calculation and visualization (make consensus member in Clustering?)
        consensus = []

        for i in range(1, max(self.clustering.max_s) + 1):
            units = self.clustering.x[np.where(self.clustering.max_s == i)[0], :]
            cluster_size = len(units[:, 0])
            units_sum = np.sum(units, axis=0)
            consensus.append(np.array([1 if x >= cluster_size / 2 else 0 for x in units_sum]))

        hc_result = linkage(pdist(consensus, metric='hamming'), method='ward')
        plt.figure(figsize=(18, 10))
        dendrogram(hc_result)
        plt.show()
