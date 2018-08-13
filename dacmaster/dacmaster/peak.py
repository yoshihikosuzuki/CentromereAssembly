import os
import sys
import random
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from interval import interval
from logzero import logger
from multiprocessing import Pool

import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go

from .io import load_unit_fasta
from .clustering import ClusteringSeqs

from BITS.utils import run_command, revcomp
from BITS.core import run_edlib

plt.style.use('ggplot')


class Peak:
    def __init__(self, index, N, unit_len, density, start_len, end_len, units):
        self.index = index   # of this peak
        self.N = N   # num of merged peaks inside this peak (1 if an independent peak)
        self.unit_len = unit_len   # unit length for each merged peak   # NOTE: list
        self.density = density   # density for each merged peak   # NOTE: list
        self.start_len = start_len   # min unit length in this peak
        self.end_len = end_len   # max unit length in this peak
        self.units = units   # longer than <start_len> and shorter than <end_len>   # NOTE: pandas dataframe

    def take_consensus_cyclic_parallel(self, args):
        """
        Return consensus sequence of the given sequences using Consed.
        First one will be the seed, and cyclic alignment is used in the mapping of other ones to it.
        This requires <out_dir> for a temporary place of the Consed input file.
        """

        read_id, path_id, seqs, tmp_dir = args

        if not os.path.isdir(tmp_dir):
            run_command(f"mkdir {tmp_dir}")
        out_fname = os.path.join(tmp_dir, f"input.seq.{os.getpid()}")   # temporary file for Consed input
        seq_writer = open(out_fname, 'w')

        # Output the first sequence as seed for consensus
        seed_seq = seqs[0]
        seq_writer.write(f"{seed_seq}\n")

        for seq in seqs[1:]:
            seq = seq * 2   # duplicated target sequence for a proxy of cyclic alignment
            alignment = run_edlib(seed_seq, seq, mode="glocal")
            #print(alignment["start"], alignment["end"], len(seq)/2)
            seq = seq[alignment["start"]:alignment["end"]]   # mapped area
            seq_writer.write(f"{seq}\n")

        seq_writer.close()
        logger.debug(f"Wrote seqs: {read_id}({path_id})")   # XXX: inserting this makes code works?

        consensus_seq = run_command(f"consed {out_fname}").replace('\n', '')
        if len(consensus_seq) == 0 or consensus_seq[0] == 'W' or consensus_seq[0] == '*':
            logger.debug(f"Strange Consed output: {read_id}({path_id})")
            return None
        else:
            logger.debug(f"Consed successed: {read_id}({path_id})")
            return (read_id, path_id, consensus_seq)

    def take_intra_consensus_parallel(self, min_n_units, n_core):
        exe_pool = Pool(n_core)

        # TODO: add a filter of minimum coverage in the read

        # TODO: divide into homogeneous TR and heterogeneous TR....

        tasks = [[read_id, path_id, list(df_path["sequence"]), "tmp"]
                 for read_id, df_read in self.units.groupby("read_id")
                 for path_id, df_path in df_read.groupby("path_id")
                 if len(df_path) >= min_n_units]

        logger.debug("tasks generated")

        self.units_consensus = {}   # intra-TR consensus units
        index = 0
        for ret in exe_pool.imap(self.take_consensus_cyclic_parallel, tasks):
            if ret is not None:
                self.units_consensus[index] = ret
                index += 1

        logger.debug("tasks completed and gathered")
        exe_pool.close()

        self.units_consensus = pd.DataFrame.from_dict(self.units_consensus,
                                                      orient="index",
                                                      columns=("read_id",
                                                               "path_id",
                                                               "sequence"))

    def construct_master_units(self, min_n_units=10, n_core=1):
        """
        From all units in a peak, construct a set of "master units".
        The main purpose of this task rather than direct construction of
        "representative units" is to fix the boundaries of the raw units in
        PacBio reads. In other words, master units have a role of "phase
        adjustment" among the units.

        However, every two units are not necessarily similar to each other.
        Since phase adjustment would become nonsense if intra-"master class"
        sequence diversity is too high, we must construct not a single master
        unit but a set of several master units in general. We expect that each
        master unit roughly corresponds to each chromosome in the genome.
        Therefore, we can say that construction of master units is most basic
        (approximately chromosome-level) characterization of the tandem repeat
        units.
        """

        # Calculate intra-TR consensus unit for each TR
        logger.debug("Start taking intra-TR consensus")
        self.take_intra_consensus_parallel(min_n_units=min_n_units, n_core=n_core)
        logger.debug("Ended intra consensus")

        # Cluster the intra-TR consensus untis
        self.clustering = ClusteringSeqs(self.units_consensus["sequence"])
        self.clustering.cluster_hierarchical(n_core=n_core)
        #self.clustering.cluster_greedy()   # too greedy way


class PeaksFinder:
    def __init__(self,
                 unit_fasta,
                 min_len=50,   # peak length must be longer than this
                 max_len=1000,   # must be shorter than this
                 band_width=5,   # parameter for KDE
                 min_density=0.001,   # threshold of peak hight
                 deviation=0.1):   # units inside of "peak_len +- deviation %" are collected as Peak.units

        if min_len < 50:
            logger.warn(f"Specified minimum unit length ({min_len} bp) is "
                        f"shorter than 50 bp, which is generally detection "
                        f"threshold of datander & datruf!")

        self.units = load_unit_fasta(unit_fasta)
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
