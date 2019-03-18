from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
from interval import interval
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
import plotly.offline as py
import plotly.graph_objs as go
from BITS.run import run_edlib
from BITS.seq import load_fasta
from BITS.utils import print_log, NoDaemonPool
import consed
from .clustering import ClusteringSeqs


@dataclass(repr=False, eq=False)
class PeaksFinder:
    """
    Class for identifying peaks in raw unit length distribution using kernel density estimation.
    """

    reads_fname: str   # DataFrame generated in run.py
    units_fname: str   # DataFrame, output of datruf
    min_len: int = 50   # only peaks longer than <min_len> and shorter than <max_len> will be found
    max_len: int = 500
    cov_th: float = 0.8   # only reads whose <cov_th> * 100 percent is covered by TR are handled
    band_width: int = 5   # param. for KDE, critical for the number of peaks detected
    min_density: float = 0.005   # threshold for peaks
    deviation: float = 0.08   # <peak_len> * (1 +- <deviation>) will be the range of each peak

    reads: pd.DataFrame = field(init=False)   # {dbid: sequence} of only TR reads
    units: pd.DataFrame = field(init=False)   # whole raw units data
    dens: np.ndarray = field(init=False)   # estimated density for each unit length
    peaks: List[Peak] = field(init=False)

    def __post_init__(self):
        if self.min_len < 50:
            logger.warn(f"You specified the minimum unit length shorter than 50 bp ({self.min_len} bp), "
                        f"which is detection limit of datander.")

        # Load all reads and all raw units
        self.reads = load_fasta(self.reads_fname)
        self.reads = {i + 1: [header, len(seq), seq]
                      for i, (header, seq) in enumerate(self.reads.items())}
        self.reads = pd.DataFrame.from_dict(self.reads,
                                            orient="index",
                                            columns=("header",
                                                     "length",
                                                     "sequence"))
        self.reads.index.names = ["dbid"]
        self.units = pd.read_csv(self.units_fname, sep='\t', index_col=0)

        # Extract TR-reads and units inside those reads
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

        original = go.Histogram(x=(self.units["length"] if entire
                                   else (self.units["length"]
                                         .pipe(lambda s: s[self.min_len < s])
                                         .pipe(lambda s: s[s < self.max_len]))),
                                xbins=dict(start=0 if entire else self.min_len,
                                           end=self.units["length"].max() if entire else self.max_len,
                                           size=1))
        layout = go.Layout(xaxis=dict(title='Unit length'),
                           yaxis=dict(title='Frequency'))
        py.iplot(go.Figure(data=[original], layout=layout))
        smoothed = go.Scatter(x=np.arange(self.min_len, self.max_len + 1),
                              y=self.dens)
        layout = go.Layout(xaxis=dict(title='Unit length'),
                           yaxis=dict(title='Density'))
        py.iplot(go.Figure(data=[smoothed], layout=layout))

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
