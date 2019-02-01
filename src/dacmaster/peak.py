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
    intvl: interval()   # of unit length of this peak

    n_sub_peaks: int = field(init=False)   # number of sub-peaks
    lens: List[int] = field(init=False)   # (sub-)peak unit lengths
    dens: List[float] = field(init=False)   # (sub-)peak densities

    def __post_init__(self, length, density):
        self.n_sub_peaks = 1
        self.lens = [length]
        self.dens = [density]

    @property
    def min_len(self):
        return int(self.intvl[0][0])

    @property
    def max_len(self):
        return int(self.intvl[0][1])

    def add_peak(self, length, density, intvl):
        self.n_sub_peaks += 1
        self.lens.append(length)
        self.dens.append(density)
        self.intvl |= intvl
        assert len(self.intvl) == 1, "Sub-peak unit length intervals must overlap"


def __take_intra_consensus(args):
    read_id, path_id, seqs = args
    cons_seq = consed.consensus([seq if i == 0
                                 else run_edlib(seqs[0],
                                                seq,
                                                "global",
                                                cyclic=True,
                                                return_seq=True).seq
                                 for i, seq in enumerate(seqs)],
                                n_iter=2)

    if cons_seq == "":
        logger.warn(f"Could not take consensus @ {read_id}({path_id})")
    else:
        logger.debug(f"Finished @ {read_id}({path_id})")

    return (read_id, path_id, len(cons_seq), cons_seq)


def _take_intra_consensus(args_list):
    return [__take_intra_consensus(args) for args in args_list]


@dataclass(repr=False, eq=False)
class Peak:
    """
    Instance variables:
      @ reads <pd.df>: all TR-reads
          Index: dbid
          [header, length, sequence]

      @ raw_units <pd.df>: unsynchronized raw units
          [read_id, path_id, start, end, length, sequence]

      @ cons_units <pd.df>: unsynchronized intra-TR consensus units
          [read_id, path_id, length, sequence]

      @ repr_units <pd.df>: synchronized, global-level representative units
          [repr_id, cluster_id, cluster_size, length, sequence]

      @ cl_master <ClusteringSeqs>: perform clustering of <raw_units> to construct master units
    """

    info: PeakInfo
    reads: pd.DataFrame
    raw_units: pd.DataFrame

    def calc_repr_units(self, min_n_units, n_core):
        """
        Main rountine.
        """

        self.take_intra_consensus(min_n_units, n_core)
        self.cluster_cons_units(n_core)
        self.generate_master_units()

    @print_log("intra-TR consensus")
    def take_intra_consensus(self, min_n_units, n_core):
        # TODO: divide into homogeneous TR and heterogeneous TR?

        tasks = [(read_id, path_id, list(df_path["sequence"]))
                 for read_id, df_read in self.raw_units.groupby("read_id")
                 for path_id, df_path in df_read.groupby("path_id")
                 if df_path.shape[0] >= min_n_units]   # filter by min. num. of units in a TR

        n_tasks = len(tasks)
        n_sub = -(-n_tasks // n_core)   # num. of tasks for each core
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
        logger.info(f"{len(self.cons_units)} out of {n_tasks} succeeded")

        self.cons_units = pd.DataFrame.from_dict(self.cons_units,
                                                 orient="index",
                                                 columns=("read_id",
                                                          "path_id",
                                                          "length",
                                                          "sequence"))

    @print_log("hierarchical clustering of consensus units")
    def cluster_cons_units(self, n_core):
        """
        Cluster the intra-TR consensus units
        """

        self.cl_master = ClusteringSeqs(self.cons_units["sequence"])
        self.cl_master.calc_dist_mat(n_core)
        self.cl_master.cluster_hierarchical()

    @print_log("master units construction")
    def generate_master_units(self,
                              redundant_threshold=0.05,
                              noisy_threshold=0.01,
                              similar_threshold=0.3):
        """
        Take consensus for each cluster of consensus units while removing noisy results.
        """

        self.master_units = self.cl_master.generate_consensus()
        logger.info(f"Original master units:\n{self.master_units}")

        # Merge too close clusters
        n_master = self.master_units.shape[0]
        for i in range(n_master - 1):
            seq_i = self.master_units["sequence"].iloc[i]
            if seq_i == "":
                continue
            for j in range(i + 1, n_master):
                seq_j = self.master_units["sequence"].iloc[j]
                if seq_j == "":
                    continue
                diff = run_edlib(seq_i,
                                 seq_j,
                                 "global",
                                 cyclic=True,
                                 rc=True,
                                 only_diff=True)
                if diff < redundant_threshold:
                    self.cl_master.merge_cluster(self.master_units["cluster_id"].iloc[i],
                                                 self.master_units["cluster_id"].iloc[j])
        self.master_units = self.cl_master.generate_consensus()
        logger.info(f"After merging close units:\n{self.master_units}")

        # Remove remaining noisy clusters
        del_row = [index for index, df in self.master_units.iterrows()
                   if df["cluster_size"] < self.cl_master.N * noisy_threshold   # too small cluster
                   or len(df["sequence"]) == 0]   # Consed failure
        self.master_units = self.master_units.drop(del_row).reset_index(drop=True)
        logger.debug(f"After removing noisy clusters:\n{self.master_units}")

        # Synchronize phase, i.e. start position
        del_row = []
        n_master = self.master_units.shape[0]
        for i in range(n_master - 1):   # TODO: simultaneously synchronize, or fix single seed
            seq_i = self.master_units["sequence"].iloc[i]
            for j in range(i + 1, n_master):
                seq_j = self.master_units["sequence"].iloc[j]
                align = run_edlib(seq_i,
                                  seq_j,
                                  "global",
                                  cyclic=True,
                                  rc=True,
                                  return_seq=True,
                                  return_seq_diff_th=similar_threshold)
                if align.seq is not None:
                    #logger.debug(f"Synchronize {i} and {j} (strand = {align['strand']})")
                    self.master_units.loc[j, "sequence"] = align.seq
        self.master_units = self.master_units.drop(del_row).reset_index(drop=True)
        logger.info(f"Final mater units:\n{self.master_units}")


@dataclass(repr=False, eq=False)
class PeaksFinder:
    """
    Class for identifying peaks in raw unit length distribution using kernel density estimation.
    """

    reads_fname: str   # DataFrame generated in run.py
    units_fname: str   # DataFrame, output of datruf
    min_len: int = 50   # only peaks longer than <min_len> and shorter than <max_len> will be found
    max_len: int = 1000
    cov_th: float = 0.8   # only reads whose <cov_th> * 100 percent is covered by TR are handled
    band_width: int = 5   # param. for KDE, critical for the number of peaks detected
    min_density: float = 0.001   # threshold for peaks
    deviation: float = 0.1   # <peak_len> * (1 +- <deviation>) will be the range of each peak

    reads: pd.DataFrame = field(init=False)   # {dbid: sequence} of only TR reads
    units: pd.DataFrame = field(init=False)   # whole raw units data
    dens: np.ndarray = field(init=False)   # estimated density for each unit length
    peaks: List[Peak] = field(init=False)

    def __post_init__(self):
        if self.min_len < 50:
            logger.warn(f"You specified the minimum unit length shorter than 50 bp ({self.min_len} bp), "
                        f"which is detection limit of datander.")

        # Extract TR-reads and units inside them
        self.reads = pd.read_csv(self.reads_fname, sep='\t', index_col=0)
        self.units = pd.read_csv(self.units_fname, sep='\t', index_col=0)
        tr_reads = set()   # read dbids to be extracted
        del_row = set()   # row indices of <self.units> to be removed
        for read_id, df in self.units.groupby("read_id"):
            if df["length"].sum() >= self.reads.loc[read_id]["length"] * self.cov_th:
                tr_reads.add(read_id)
            else:
                del_row.update(df.index)

        n_all_reads, n_all_units = self.reads.shape[0], self.units.shape[0]
        self.reads = self.reads.loc[sorted(tr_reads)]
        self.units = self.units.drop(del_row).reset_index(drop=True)
        logger.info(f"{self.reads.shape[0]} out of {n_all_reads} all reads and "
                    f"{self.units.shape[0]} out of {n_all_units} all units were loaded")

        # Smooth the unit length distribution by KDE
        self.smooth_dist()

    def smooth_dist(self):
        """
        Run KDE and calculate density for each unit length value.
        """

        self.dens = np.exp(KernelDensity(kernel='gaussian',
                                         bandwidth=self.band_width)
                           .fit(np.array(self.units["length"]
                                         .pipe(lambda s: s[self.min_len < s])
                                         .pipe(lambda s: s[s < self.max_len]))
                                .reshape(-1, 1))
                           .score_samples(np.arange(self.min_len, self.max_len + 1)
                                          .reshape(-1, 1)))
        assert self.dens.shape[0] == (self.max_len - self.min_len + 1), "Inconsistent data range"

    def plot_dist(self, entire=False):
        """
        Show both original unit length distribution and smoothed distribution.
        """

        original = go.Histogram(x=(self.units["length"]
                                   .pipe(lambda s: s[self.min_len < s])
                                   .pipe(lambda s: s[s < self.max_len])),
                                xbins=dict(start=0 if entire else self.min_len,
                                           end=self.units["length"].max() if entire else self.max_len,
                                           size=1))
        smoothed = go.Scatter(x=np.arange(self.min_len, self.max_len + 1),
                              y=self.dens,
                              yaxis='y2')
        layout = go.Layout(xaxis=dict(title='Unit length'),
                           yaxis=dict(title='Frequency'),
                           yaxis2=dict(title='Density', side='right'))
        py.iplot(go.Figure(data=[original, smoothed], layout=layout))

    @print_log("peak detection")
    def run(self):
        """
        Detect peaks in the unit length distribution.
        Adjascent peaks close to each other are merged into single peak.
        """

        # Naive sweep
        peak_infos = []
        for i in range(1, self.max_len - self.min_len):
            if not ((self.dens[i] > self.dens[i - 1]) and
                    (self.dens[i] > self.dens[i + 1]) and
                    (self.dens[i] >= self.min_density)):
                continue

            l, d = self.min_len + i, self.dens[i]   # peak unit length and its density
            intvl = interval[-(- l * (1. - self.deviation) // 1),
                             int(l * (1. + self.deviation))]

            # Add a peak
            if len(peak_infos) == 0 or peak_infos[-1].intvl & intvl == interval():
                logger.info(f"New peak detected: {l} bp (density = {d:.5f})")
                peak_infos.append(PeakInfo(l, d, intvl))
            else:
                # Merge close peaks
                logger.info(f"Sub-peak merged: {l} bp (density = {d:.5f})")
                peak_infos[-1].add_peak(l, d, intvl)

        self.peaks = []
        for peak_info in peak_infos:
            # Extract units and reads belonging to the peak
            raw_units = (self.units[self.units["length"] >= peak_info.min_len]
                         .pipe(lambda df: df[df["length"] <= peak_info.max_len]))
            reads = self.reads.loc[set(raw_units["read_id"])]
            self.peaks.append(Peak(peak_info, reads, raw_units))
