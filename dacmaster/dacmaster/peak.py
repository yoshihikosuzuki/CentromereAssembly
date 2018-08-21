import os
import copy
import numpy as np
import pandas as pd
import networkx as nx
from interval import interval
from logzero import logger
from multiprocessing import Pool

import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go

from .clustering import ClusteringSeqs

#from BITS.utils import run_command
from BITS.seq import revcomp
from BITS.run import run_edlib
import consed

plt.style.use('ggplot')


def __take_intra_consensus(args):
    read_id, path_id, seqs = args
    cons_seq = consed.consensus([seq if i == 0
                                 else run_edlib(seqs[0],
                                                seq,
                                                mode="glocal",
                                                cyclic=True,
                                                return_seq=True)["seq"]
                                 for i, seq in enumerate(seqs)])

    if cons_seq == "":
        logger.warn(f"Could not take consensus @ {read_id}({path_id})")
    else:
        logger.debug(f"Finished @ {read_id}({path_id})")

    return (read_id, path_id, cons_seq)


def _take_intra_consensus(args_list):
    return [__take_intra_consensus(args) for args in args_list]


class Peak:
    """
    Description of the instance variables:
      <peak_id>      : int              : peaks[peak_id] should be itself in run.py
      <N>            : int              : num of merged peaks inside this peak (1 if an isolated single peak)
      <unit_len>     : list(int)        : unit length for each merged peak
      <density>      : list(float)      : density for each merged peak
      <min_len>      : int              : min unit length in this peak
      <max_len>      : int              : max unit length in this peak

      <reads>        : pd.df            : # TODO: let this include encodings?
        [read_id, sequence]

      <raw_units>    : pd.df            : unsynchronized raw units between <min_len> and <max_len> bp
        [read_id, path_id, start, end, length]
      <cons_units>   : pd.df            : unsynchronized intra-TR consensus units
        [read_id, path_id, length, sequence]
      <master_units> : pd.df            : synchronized, global-level representative units
        [master_id, cluster_id, cluster_size, n_bad_align, length, sequence]
      <repr_units>   : pd.df            : synchronized, local-level representative units
        [repr_id, cluster_id, cluster_size, length, sequence]
      <encodings>    : pd.df            : synchronized raw units   # TODO: separate into master_encodings and repr_encodings?

      <cl_master>    : ClusteringSeqs   : perform clustering of <raw_units> to construct master units
      <cl_repr>      : ClusteringSeqs   : 
      <cl_unit>      : ClusteringVarMat :
    """

    # TODO:
    # func to load <reads>
    # remove unit_id and sequence columns in <raw_units>

    def __init__(self, peak_id, N, unit_len, density, min_len, max_len, raw_units):
        self.peak_id = peak_id
        self.N = N
        self.unit_len = unit_len
        self.density = density
        self.min_len = min_len
        self.max_len = max_len
        self.raw_units = raw_units

    def take_intra_consensus(self, min_n_units, n_core):
        # TODO: divide into homogeneous TR and heterogeneous TR?

        tasks = [(read_id, path_id, list(df_path["sequence"]))
                 for read_id, df_read in self.raw_units.groupby("read_id")
                 for path_id, df_path in df_read.groupby("path_id")
                 if df_path.shape[0] >= min_n_units]   # filter by min. num. of units in a TR   # TODO: add a filter of minimum coverage in the read

        n_sub = -(-len(tasks) // n_core)   # num. of tasks for each core
        tasks_sub = [tasks[i * n_sub:(i + 1) * n_sub - 1] for i in range(n_core)]

        self.cons_units = {}
        index = 0
        logger.debug(f"Scattering tasks with {n_core} cores")
        exe_pool = Pool(n_core)
        for ret in exe_pool.imap(_take_intra_consensus, tasks_sub):
            logger.debug(f"Received")
            for r in ret:
                self.cons_units[index] = r
                index += 1
        exe_pool.close()
        logger.debug("Finished all tasks")

        self.cons_units = pd.DataFrame.from_dict(self.cons_units,
                                                 orient="index",
                                                 columns=("read_id",
                                                          "path_id",
                                                          "sequence"))

    def adjust_two_units(self, seed, seq):   # allow only forward-forward
        return run_edlib(seed, seq, mode="glocal", cyclic=True, return_seq=True)["seq"]

    def filter_master_units(self, redundant_threshold=0.05, similar_threshold=0.2):
        """
        Filter "noisy" master units from the cluster consensus sequences.
        """

        self.master_units = copy.deepcopy(self.master_original)

        del_row = []
        for index, df in self.master_original.iterrows():
            cluster_id, cluster_size, skip_count, seq = df
            if cluster_size < self.cl_master.N * 0.01:   # too small cluster
                del_row.append(index)
            elif len(seq) == 0 or seq[0] in ['W', '*']:   # strange consensus seq
                del_row.append(index)
        self.master_units = self.master_units.drop(del_row)
        self.master_units = self.master_units.reset_index(drop=True)

        logger.debug(f"After removing noisy clusters:\n{self.master_units}")
        
        # due to redundancy
        n_master = self.master_units.shape[0]
        dm = np.zeros((n_master, n_master), dtype='float32')
        strand = np.zeros((n_master, n_master), dtype='int')    # 0 for forward, 1 for revcomp
        del_row = []
        for i in range(n_master - 1):
            for j in range(i + 1, n_master):
                align = run_edlib(self.master_units["sequence"][i],
                                  self.master_units["sequence"][j],
                                  mode="glocal",
                                  revcomp=True,
                                  cyclic=True)
                dm[i, j] = dm[j, i] = align["diff"]
                strand[i, j] = strand[j, i] = align["strand"]
                if dm[i, j] < redundant_threshold:
                    if self.master_units["cluster_size"][i] >= self.master_units["cluster_size"][j]:
                        del_row.append(j)
                    else:
                        del_row.append(i)
        logger.debug(f"\ndist mat:\n{dm}\nstrand:\n{strand}")
        self.master_units = self.master_units.drop(del_row)
        self.master_units = self.master_units.reset_index(drop=True)
        np.delete(dm, del_row, axis=0)
        np.delete(dm, del_row, axis=1)
        np.delete(strand, del_row, axis=0)
        np.delete(strand, del_row, axis=1)
        logger.debug(f"After removing redundancy:\n{self.master_units}")
        logger.debug(f"\ndist mat:\n{dm}\nstrand:\n{strand}")

        # for similar seqs, adjust start positions and strands
        for i in range(dm.shape[0]):
            dm[i, i] = np.inf
        g = nx.Graph()
        while np.min(dm) <= similar_threshold:
            i, j = np.unravel_index(dm.argmin(), dm.shape)   # index of min diff
            logger.debug(f"argmin @ {i}, {j}")
            g.add_node(i)
            g.add_node(j)
            g.add_edge(i, j)
            dm[i, j] = dm[j, i] = np.inf
        for c in nx.connected_components(g):
            nodes = sorted([node for node in c])
            logger.debug(f"cc: {nodes}")
            seed = nodes[0]
            revs = [node for node in nodes if strand[seed, node] == 1]
            logger.debug(f"rev: {revs}")
            for rev in revs:   # strand
                self.master_units.loc[rev, "sequence"] = revcomp(self.master_units.loc[rev, "sequence"])
                strand[rev, :] = 1 - strand[rev, :]
                strand[:, rev] = 1 - strand[:, rev]
            for node in nodes[1:]:   # revcomp
                self.master_units.loc[node, "sequence"] = self.adjust_two_units(self.master_units.loc[seed, "sequence"],
                                                                                self.master_units.loc[node, "sequence"])   # TODO: simply running and replacing to run_edlib(return_seq=True, return_seq_diff_th=similar_threshold)["seq"] for each row is enough?
        logger.debug(f"\n{self.master_units}\ndist mat:\n{dm}\nstrand:\n{strand}")

    def construct_master_units(self, min_n_units, n_core):
        """
        From all units in a peak, construct a set of "master units".
        The main purpose here is to adjust and fix the boundaries of the raw
        units in the reads. In other words, master units play a role of "phase
        adjustment" among the units. As side effect, these master units play
        another role of large scale-level (almost chromosome-level) segregation
        of the units.
        """

        # Calculate intra-TR consensus unit for each TR
        if not hasattr(self, "cons_units"):
            logger.info("Starting taking intra-TR consensus")
            self.take_intra_consensus(min_n_units, n_core)
            logger.info("Finished intra-TR consensus")

        # Cluster the intra-TR consensus units
        if not hasattr(self, "cl_master"):
            self.cl_master = ClusteringSeqs(self.cons_units["sequence"])
        logger.info("Starting hierarchical clustering")
        self.cl_master.cluster_hierarchical(n_core=n_core)
        #self.cl_master.cluster_greedy()   # too greedy way
        logger.info("Finished hierarchical clustering")

        # Generate master units
        self.master_original = self.cl_master.generate_consensus()   # NOTE: dataframe
        logger.debug(f"\n{self.master_original}")
        # remove strange sequences and redundancy
        # after that, start-position adjustment and strand adjustment of close seqs
        self.filter_master_units()
        logger.info("Finished cluster consensus")


class PeaksFinder:
    def __init__(self,
                 units_fname,
                 min_len=50,   # peak length must be longer than this
                 max_len=1000,   # must be shorter than this
                 band_width=5,   # parameter for KDE
                 min_density=0.001,   # threshold of peak hight
                 deviation=0.1):   # units inside of "peak_len +- deviation %" are collected as Peak.units

        if min_len < 50:
            logger.warn(f"Specified minimum unit length ({min_len} bp) is "
                        f"shorter than 50 bp, which is generally detection "
                        f"threshold of datander & datruf!")

        self.units = pd.read_table(units_fname, index_col=0)
        self.min_len = min_len
        self.max_len = max_len
        self.band_width = band_width
        self.min_density = min_density
        self.deviation = deviation

    def run(self):
        self.smooth_unit_len_dist()
        self.detect_peaks()

    def smooth_unit_len_dist(self):
        """
        Calculate unit length distribution smoothed by kernel density estimation.
        """

        from sklearn.neighbors import KernelDensity

        self.unit_lens = np.array(self.units["length"])

        # all unit lengths within the specified interval
        self.ul = self.unit_lens[(self.min_len < self.unit_lens)
                                 & (self.unit_lens < self.max_len)]

        # [min_len, min_len + 1, ..., max_len]
        self.ls = np.linspace(self.min_len,
                              self.max_len,
                              self.max_len - self.min_len + 1,
                              dtype=int)

        KDE = KernelDensity(kernel='gaussian',
                            bandwidth=self.band_width)
        self.kde = KDE.fit(self.ul.reshape(-1, 1))

        # estimated density at each unit length
        self.dens = np.exp(self.kde.score_samples(self.ls.reshape(-1, 1)))

    def detect_peaks(self):
        """
        Detect peaks in the unit length distribution with a simple sweep line.
        Adjascent peaks close to each other are merged.
        """

        self.peak_intervals = interval()
        self.peak_info = []

        # First detect peaks
        prev_n_peaks = 0
        for i in range(1, len(self.dens) - 1):
            if ((self.dens[i] >= self.min_density) and
                    (self.dens[i - 1] < self.dens[i]) and
                    (self.dens[i] > self.dens[i + 1])):

                # collect units, allowing the deviation at both sides
                self.peak_intervals |= interval[-(- self.ls[i] * (1. - self.deviation) // 1),
                                                int(self.ls[i] * (1. + self.deviation))]
                # also we keep information of unit length and density
                peak_info = (self.ls[i], self.dens[i])

                logger.info(f"Peak detected: length = {self.ls[i]} bp, density = {self.dens[i]}")

                if prev_n_peaks == len(self.peak_intervals):
                    # the peak detected in this loop has been merged to the previous one
                    self.peak_info[-1].append(peak_info)
                    logger.info("Merged to the previous peak.")
                else:
                    # new peak interval is independent
                    self.peak_info.append([peak_info])
                    prev_n_peaks += 1

        # For each peak interval, create Peak class instance
        self.peaks = []
        for i, intvl in enumerate(self.peak_intervals.components):
            peak_info = self.peak_info[i]
            N = len(peak_info)   # num of peaks merged
            length, density = list(zip(*peak_info))   # list for each merged peak
            start, end = intvl[0]   # min peak_len - deviation, max peak_len + deviation
            units = (self.units[self.units["length"] >= start]   # belonging to this peak
                     .pipe(lambda df: df[df["length"] <= end]))

            self.peaks.append(Peak(i, N, length, density, start, end, units))

    def plot_unit_len_dist(self):
        # Raw unit length distribution
        data = [go.Histogram(x=self.unit_lens,
                             xbins=dict(start=self.min_len,
                                        end=self.max_len,
                                        size=1))]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Frequency"))
        py.iplot(go.Figure(data=data, layout=layout))

        # Smoothed distribution
        data = [go.Scatter(x=self.ls, y=self.dens)]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Density"))
        py.iplot(go.Figure(data=data, layout=layout))
