import pickle
from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
from interval import interval
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
from BITS.run import run_edlib
from BITS.utils import print_log, NoDaemonPool
import consed
from .clustering import ClusteringSeqs

plt.style.use('ggplot')


@dataclass(repr=False, eq=False)
class PeakInfo:
    """
    Metadata for Peak class. Also used in the peak detection.
    """

    length: InitVar[int]
    density: InitVar[float]
    intvl: interval()   # units whose lengths are within this interval belong to this peak

    lens: List[int] = field(init=False)   # unit lengths of sub-peaks
    dens: List[float] = field(init=False)   # density of sub-peaks

    def __post_init__(self, length, density):
        self.lens = [length]
        self.dens = [density]

    @property
    def N(self):
        """
        Get the number of sub-peaks inside of this peak.
        """
        assert len(self.lens) == len(self.dens), "Inconsistent unit length range"
        return len(self.lens)

    @property
    def min_len(self):
        """
        Get the minimum unit length in this peak.
        """
        return self.intvl[0][0]

    @property
    def max_len(self):
        """
        Get the maximum unit length in this peak.
        """
        return self.intvl[0][1]

    def add_peak(self, length, density, intvl):
        self.lens.append(length)
        self.dens.append(density)
        self.intvl |= intvl
        assert len(self.intvl) == 1, "Split peak intervals"

    def overlaps_to(self, intvl):
        return True if self.intvl & intvl != interval() else False


def save_peaks(peaks, pkl_fname="peaks.pkl"):
    with open(pkl_fname, 'wb') as f:
        pickle.dump(peaks, f)


def load_peaks(pkl_fname="peaks.pkl"):
    with open(pkl_fname, 'rb') as f:
        peaks = [peak_from_old_data(p) for p in pickle.load(f)]
        logger.info(f"{len(peaks)} peaks were loaded")
        return peaks


def peak_from_old_data(p):
    """
    Reconstruct Peak instance and its instance variables with the latest instance methods.
    """

    peak = Peak(None, None)

    # instance variables of Peak class
    for attr in ("info", "raw_untis", "cons_units", "master_units", "master_original"):
        if hasattr(p, attr):
            setattr(peak, attr, getattr(p, attr))

    # ClusteringSeqs instance
    if hasattr(peak, "cons_units"):
        peak.cl_master = ClusteringSeqs(peak.cons_units["sequence"])

        # instance variables of ClusteringSeqs class
        for attr in ("hc_result", "hc_result_precomputed", "assignment", "dist_matrix", "hc_input", "coord"):
            if hasattr(p.cl_master, attr):
                setattr(peak.cl_master, attr, getattr(p.cl_master, attr))

    return peak


def __take_intra_consensus(args):
    read_id, path_id, seqs = args
    cons_seq = consed.consensus([seq if i == 0
                                 else run_edlib(seqs[0],
                                                seq,
                                                mode="glocal",
                                                cyclic=True,
                                                return_seq=True)["seq"]
                                 for i, seq in enumerate(seqs)],
                                n_iter=2)

    if cons_seq == "":
        logger.warn(f"Could not take consensus @ {read_id}({path_id})")
    else:
        logger.debug(f"Finished @ {read_id}({path_id})")

    return (read_id, path_id, cons_seq)


def _take_intra_consensus(args_list):
    return [__take_intra_consensus(args) for args in args_list]


@dataclass(repr=False, eq=False)
class Peak:
    """
    Instance variables:
      @ reads <pd.df>: # TODO: let this include encodings?
          [read_id, sequence]

      @ raw_units <pd.df>: unsynchronized raw units
          [read_id, path_id, start, end, length]

      @ cons_units <pd.df>: unsynchronized intra-TR consensus units
          [read_id, path_id, length, sequence]

      @ master_units <pd.df>: synchronized, global-level representative units
          [master_id, cluster_id, cluster_size, n_bad_align, length, sequence]

      @ repr_units <pd.df>: synchronized, local-level representative units
          [repr_id, cluster_id, cluster_size, length, sequence]

      @ encodings <pd.df>: by synchronized raw units

      @ cl_master <ClusteringSeqs>: perform clustering of <raw_units> to construct master units

      @ cl_repr <ClusteringSeqs>: 

      @ cl_unit <ClusteringVarMat>:
    """

    info: PeakInfo
    raw_units: pd.DataFrame

    @print_log("intra-TR consensus")
    def take_intra_consensus(self, min_n_units, n_core):
        # TODO: divide into homogeneous TR and heterogeneous TR?

        tasks = [(read_id, path_id, list(df_path["sequence"]))
                 for read_id, df_read in self.raw_units.groupby("read_id")
                 for path_id, df_path in df_read.groupby("path_id")
                 if df_path.shape[0] >= min_n_units]   # filter by min. num. of units in a TR

        n_sub = -(-len(tasks) // n_core)   # num. of tasks for each core
        tasks_sub = [tasks[i * n_sub:(i + 1) * n_sub - 1] for i in range(n_core)]

        self.cons_units = {}
        index = 0
        logger.debug(f"Scattering tasks with {n_core} cores")
        exe_pool = NoDaemonPool(n_core)
        for ret in exe_pool.map(_take_intra_consensus, tasks_sub):
            for r in ret:
                if r[2] != "":
                    self.cons_units[index] = r
                    index += 1
        exe_pool.close()
        exe_pool.join()

        self.cons_units = pd.DataFrame.from_dict(self.cons_units,
                                                 orient="index",
                                                 columns=("read_id",
                                                          "path_id",
                                                          "sequence"))

    @print_log("hierarchical clustering of consensus units")
    def cluster_cons_units(self, n_core):
        # Cluster the intra-TR consensus units
        if not hasattr(self, "cl_master"):
            self.cl_master = ClusteringSeqs(self.cons_units["sequence"])
        self.cl_master.cluster_hierarchical(n_core=n_core)

    @print_log("master units construction")
    def generate_master_units(self, redundant_threshold=0.05, similar_threshold=0.3):
        """
        Take consensus for each cluster of consensus units while removing noisy results.
        """

        self.master_units = self.cl_master.generate_consensus()
        logger.info(f"Original master units:\n{self.master_units}")

        # Remove noisy clusters
        del_row = [index for index, df in self.master_units.iterrows()
                   if df["cluster_size"] < self.cl_master.N * 0.01    # too small cluster
                   or len(df["sequence"]) == 0]   # Consed failure
        self.master_units = self.master_units.drop(del_row).reset_index(drop=True)
        logger.debug(f"After removing noisy clusters:\n{self.master_units}")

        # Remove redundant clusters   # TODO: recompute consensus sequence after merging similar clusters
        n_master = self.master_units.shape[0]
        del_row = [i if self.master_units["cluster_size"][i] < self.master_units["cluster_size"][j]
                   else j
                   for i in range(n_master - 1)
                   for j in range(i + 1, n_master)
                   if run_edlib(self.master_units["sequence"][i],
                                self.master_units["sequence"][j],
                                mode="glocal",
                                revcomp=True,
                                cyclic=True)["diff"] < redundant_threshold]
        self.master_units = self.master_units.drop(del_row).reset_index(drop=True)
        logger.debug(f"After removing redundancy:\n{self.master_units}")

        # Adjust start positions and strands of similar master units
        n_master = self.master_units.shape[0]
        for i in range(n_master - 1):
            for j in range(i + 1, n_master):
                align = run_edlib(self.master_units.loc[i, "sequence"],
                                  self.master_units.loc[j, "sequence"],
                                  mode="glocal",
                                  cyclic=True,
                                  revcomp=True,
                                  return_seq=True,
                                  return_seq_diff_th=similar_threshold)
                if align["seq"] is not None:
                    logger.debug(f"Synchronized {i} and {j} (strand = {align['strand']})")
                    self.master_units.loc[j, "sequence"] = align["seq"]
        logger.info(f"Final mater units:\n{self.master_units}")


@dataclass(repr=False, eq=False)
class PeaksFinder:
    """
    Class for identifying peaks in raw unit length distribution using kernel density estimation.
    """

    units_fname: InitVar[str] = "datruf_units"
    min_len: int = 50   # only peaks longer than <min_len> and shorter than <max_len> will be found
    max_len: int = 1000
    band_width: int = 5   # param. for KDE
    min_density: float = 0.001   # threshold for peaks
    deviation: float = 0.1   # <peak_len> * (1 +- <deviation>) will be the range of each peak

    units: pd.DataFrame = field(init=False)   # whole raw units data
    ulens: np.ndarray = field(init=False)   # unit lengths within [min_len, max_len]
    x: np.ndarray = field(init=False)   # [min_len, min_len + 1, ..., max_len]
    dens: np.ndarray = field(init=False)   # estimated density for each unit length
    peaks: List[Peak] = field(init=False)

    def __post_init__(self, units_fname):
        if self.min_len < 50:
            logger.warn(f"Specified minimum unit length ({self.min_len} bp) is shorter than 50 bp, "
                        f"which is typical detection limit of datander & datruf!")

        self.units = pd.read_table(units_fname, index_col=0)   # TODO: exclude "unit_id" and "sequence"
        self.ulens = np.array(self.units["length"]
                              .pipe(lambda s: s[self.min_len < s])
                              .pipe(lambda s: s[s < self.max_len]))
        self.x = np.arange(self.min_len, self.max_len + 1)

        self.smooth_dist()

    def smooth_dist(self):
        """
        Run KDE and calculate density for each unit length value.
        """

        self.dens = np.exp(KernelDensity(kernel='gaussian',
                                         bandwidth=self.band_width)
                           .fit(self.ulens.reshape(-1, 1))
                           .score_samples(self.x.reshape(-1, 1)))

    def plot_dist(self, entire=False):
        """
        Show both original unit length distribution and smoothed distribution.
        """

        # Original distribution
        start, end = ((self.min_len, self.max_len) if not entire
                      else (np.min(self.ulens), np.max(self.ulens)))
        data = [go.Histogram(x=self.ulens,
                             xbins=dict(start=start, end=end, size=1))]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Frequency"))
        py.iplot(go.Figure(data=data, layout=layout))

        # Smoothed distribution
        data = [go.Scatter(x=self.x, y=self.dens)]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Density"))
        py.iplot(go.Figure(data=data, layout=layout))

    @print_log("peak detection")
    def detect_peaks(self):
        """
        Detect peaks in the unit length distribution.
        Adjascent peaks close to each other are merged into single peak.
        """

        assert self.x.shape[0] == self.dens.shape[0], "Inconsistent variable lengths"

        # NOTE: not unit length but index on self.x and self.dens
        peak_index = [i for i in range(1, self.dens.shape[0] - 1)
                      if ((self.dens[i] > self.dens[i - 1]) and
                          (self.dens[i] > self.dens[i + 1]) and
                          (self.dens[i] >= self.min_density))]

        # Merge close peaks and prepare metadata for each peak
        peak_infos = []
        for i in peak_index:
            intvl = interval[-(- self.x[i] * (1. - self.deviation) // 1),
                             int(self.x[i] * (1. + self.deviation))]

            if len(peak_infos) == 0 or not peak_infos[-1].overlaps_to(intvl):
                logger.info(f"New peak detected: {self.x[i]} bp (density = {self.dens[i]:.5f})")
                peak_infos.append(PeakInfo(self.x[i], self.dens[i], intvl))
            else:
                logger.info(f"Sub-peak merged: {self.x[i]} bp (density = {self.dens[i]:.5f})")
                peak_infos[-1].add_peak(self.x[i], self.dens[i], intvl)

        return [Peak(peak_info,
                     (self.units[self.units["length"] >= peak_info.min_len]
                      .pipe(lambda df: df[df["length"] <= peak_info.max_len])))
                for peak_info in peak_infos]
