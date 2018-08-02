import os
import pickle
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from interval import interval
from logzero import logger

from .dacmaster_io import load_unit_fasta

from BITS.utils import run_command, revcomp
from BITS.core import run_edlib


def take_consensus_cyclic(seqs, out_dir):
    """
    Return consensus sequence of the given sequences using Consed.
    First one will be the seed, and cyclic alignment is used in the mapping of other ones to it.
    This requires <out_dir> for a temporary place of the Consed input file.
    """

    if not os.path.isdir(out_dir):
        run_command(f"mkdir {out_dir}")
    out_fname = os.path.join(out_dir, "input.seq")   # temporary file for Consed input
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
    return run_command(f"consed {out_fname}").replace('\n', '')


class Peak:
    def __init__(self, root_dir, fname_prefix, index, N, unit_len, density, start_len, end_len, units):
        self.root_dir = root_dir
        self.fname_prefix = fname_prefix
        self.index = index   # "<root_dir>/<fname_prefix>.<index>.fasta" is the result file
        self.N = N   # num of merged peaks inside this peak (1 if an independent peak)
        self.unit_len = unit_len   # unit length for each merged peak   # NOTE: list
        self.density = density   # density for each merged peak   # NOTE: list
        self.start_len = start_len   # min unit length in this peak
        self.end_len = end_len   # max unit length in this peak
        self.units = units   # longer than <start_len> and shorter than <end_len>   # NOTE: pandas dataframe

        # pickle of this class instance itself
        self.pkl_fname = os.path.join(self.root_dir,
                                      f"{self.fname_prefix}.{self.index}.pkl")

    def write(self):
        with open(self.pkl_fname, 'wb') as f:
            pickle.dump(self, f)

    def filter_units_by_single_alignment_coverage(self, threshold=0.9):
        """
        Filter units in each read where at least 100 * <threshold> % of the read
        is not covered by a single TR alignment.

        This is used to extract pure and relatively clean centromeric monomers.
        """

        # TODO: this should be done in datruf!
        
    def take_intra_consensus(self, n_units_threshold=5):
        """
        For each tandem repeat that has at least <n_units_threshold> units,
        this calculates intra-TR consensus sequence of the units.
        """

        self.units_consensus = {}   # intra-TR consensus units
        index = 0
        for read_id, df_read in self.units.groupby("read_id"):   # for each read
            for path_id, df_path in df_read.groupby("path_id"):   # for each TR
                if len(df_path) < n_units_threshold:   # only small number of units inside the TR
                    continue
                consensus_seq = take_consensus_cyclic(list(df_path["sequence"]), "tmp")
                #print(consensus_seq)
                if len(consensus_seq) == 0 or consensus_seq[0] == 'W' or consensus_seq[0] == '*':
                    # XXX: TODO: this is just a workaround for the output of consed "Warning: tile overlaps did not align well" or "** EXISTING". Fix it.
                    continue
                self.units_consensus[index] = [read_id,
                                               path_id,
                                               df_path["start"].iloc[0],   # TODO: maybe no need
                                               df_path["end"].iloc[-1],
                                               len(consensus_seq),
                                               consensus_seq]
                index += 1

        self.units_consensus = pd.DataFrame.from_dict(self.units_consensus,
                                                      orient="index",
                                                      columns=("read_id",
                                                               "path_id",
                                                               "start",
                                                               "end",
                                                               "length",
                                                               "sequence"))

    def find_representatives(self):
        """
        From all units in a peak, find representative monomers that do not align well to each other.
        This task includes start-position adjustment and phase-adjustment as well.

        In this step, determination of seed units is essential.
        """

        # Calculate consensus of the units for each TR
        self.take_intra_consensus()

        # Clustering of the consensus untis
        # Since what we are doing for now is just determine rough representative monomers,
        # this task is rather just to eliminate redundant consensus units.
        

class Runner:
    def __init__(self,
                 unit_fasta,
                 peaks_dir,
                 peak_fname_prefix,
                 min_len=50,   # peak length must be longer than this
                 max_len=1000,   # must be shorter than this
                 band_width=5,   # parameter for KDE
                 min_density=0.001,   # threshold of peak hight
                 deviation=0.1):   # units inside of "peak_len +- deviation %" are collected as Peak.units

        self.units = load_unit_fasta(unit_fasta)
        self.peaks_dir = peaks_dir
        self.peak_fname_prefix = peak_fname_prefix
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

        self.ul = self.unit_lens[(self.min_len < self.unit_lens) & (self.unit_lens < self.max_len)]
        self.ls = np.linspace(self.min_len, self.max_len, self.max_len - self.min_len + 1, dtype=int)

        self.kde = KernelDensity(kernel='gaussian', bandwidth=self.band_width).fit(self.ul.reshape(-1, 1))
        self.dens = np.exp(self.kde.score_samples(self.ls.reshape(-1, 1)))

    def detect_peaks(self):
        """
        Detect peaks in the unit length distribution by simple sweep line.
        Adjascent peaks close to each other will be automatically merged.
        """

        self.peak_intervals = interval()
        self.peak_info = []
        prev_n_peaks = 0
        for i in range(1, len(self.dens) - 1):
            if ((self.dens[i] >= self.min_density) and
                (self.dens[i - 1] < self.dens[i]) and
                (self.dens[i] > self.dens[i + 1])):

                self.peak_intervals |= interval[-(- self.ls[i] * (1. - self.deviation) // 1),
                                           int(self.ls[i] * (1. + self.deviation))]
                peak_info = (self.ls[i], self.dens[i])

                logger.info(f"Peak detected: length = {self.ls[i]} bp, density = {self.dens[i]}")
                if prev_n_peaks == len(self.peak_intervals):
                    self.peak_info[-1].append(peak_info)
                    logger.info("Merged to the previous peak.")
                else:   # new peak interval is independent
                    self.peak_info.append([peak_info])
                    prev_n_peaks += 1
    
        # For each peak interval, create Peak class instance
        self.peaks = []
        for i, intvl in enumerate(self.peak_intervals.components):
            peak_info = self.peak_info[i]
            peak_n = len(peak_info)   # num of peaks merged
            peak_len, peak_dens = list(zip(*peak_info))   # list for each merged peak
            peak_start, peak_end = intvl[0]   # min peak_len - deviation, max peak_len + deviation
            peak_units = (self.units[self.units["length"] >= peak_start]
                          .pipe(lambda df: df[df["length"] <= peak_end]))
            self.peaks.append(Peak(self.peaks_dir,
                                   self.peak_fname_prefix,
                                   i,
                                   peak_n,
                                   peak_len,
                                   peak_dens,
                                   peak_start,
                                   peak_end,
                                   peak_units))

























"""
    def align_start_positions(self):
        # For each peak
        for peak_id in range(1, self.n_peaks + 1):
            # NOTE: insert any filter on peak_id here
            if peak_id != 3:
                continue

            self.align_start_positions_peak(peak_id)

    def align_start_positions_peak(self, peak_id):
        # For each unit in the peak
        for unit_id in range(len(self.peak_units[peak_id])):
            # NOTE: insert any filter on unit_id here
            if unit_id != 1:
                continue

            self.align_start_positions_peak_unit(peak_id, unit_id)

    def align_start_positions_peak_unit(self, peak_id, unit_id, max_diff=0.28):
        out_fname_prefix = ("%s/%s.%d.aligned.%d"
                            % (self.peaks_dir,
                               self.peak_fname_prefix,
                               peak_id,
                               unit_id))

        seq_writer = open(out_fname_prefix + ".seq", 'w')   # input of Consed
        fasta_writer = open(out_fname_prefix + ".fasta", 'w')   # for split_dis
        sam_writer = open(out_fname_prefix + ".sam", 'w')   # for split_dis

        # used for reverse complement
        RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

        # a (global) -> b (local)
        def align_edlib(a, b):
            ret = edlib.align(a, b, mode="HW", task="path")

            cigar = ret["cigar"]
            # Swap I and D in cigar to satisfy the definition in SAM format
            # NOTE: This is because: in each mapping the dup is reference,
            # but indeed the unit is overall the seed (= reference) in SAM
            reverse_cigar = cigar.replace("I", "?").replace("D", "I").replace("?", "D")

            alignment_len = sum(list(map(int, re.findall(r'([0-9]+)', cigar))))
            diff = ret["editDistance"] / alignment_len   # 0 ~ 1
            start, end = ret["locations"][0]
            end += 1

            return {"cigar": reverse_cigar,
                    "diff": diff,
                    "start": start,
                    "end": end}

        unit_header, unit_seq = self.peak_units[peak_id][unit_id]

        # Output the seed unit into fasta and sam
        seq_writer.write("%s\n" % unit_seq)
        fasta_writer.write(">%s\n%s\n" % (unit_header, unit_seq))
        sam_writer.write("@SQ\tSN:%s\tLN:%d\n" % (unit_header, len(unit_seq)))

        # TODO: store strand information for assembly

        for unit in self.peak_units[peak_id]:
            dup_header = unit[0]

            # The seed unit should be already written in the first line
            if dup_header == unit_header:
                continue

            dup_seq = unit[1] * 2

            f_map = align_edlib(unit_seq, dup_seq)

            rc_seq = "".join([RC_MAP[c] for c in dup_seq[::-1]])
            rc_map = align_edlib(unit_seq, rc_seq)

            # Filter bad alignments because Consed failed if a unit cannot be aligned to the seed
            if min(f_map["diff"], rc_map["diff"]) > max_diff:
                continue

            if f_map["diff"] <= rc_map["diff"]:
                out_map = f_map
                out_seq = dup_seq
            else:
                out_map = rc_map
                out_seq = rc_seq

            out_header = "%s/%d_%d" % (dup_header.split("0_")[0],
                                       out_map["start"],
                                       out_map["end"])
            out_seq = out_seq[out_map["start"]:out_map["end"]]

            seq_writer.write("%s\n" % out_seq)
            fasta_writer.write(">%s\n%s\n" % (out_header, out_seq))
            sam_writer.write("\t".join(list(map(str, [out_header,   # QNAME
                                                      0,   # (alignment status) FLAG
                                                      unit_header,   # RNAME
                                                      1,   # (start) POS (in reference; 1-index)
                                                      30,   # MAPQ
                                                      out_map["cigar"],   # CIGAR (R->Q)
                                                      "*",   # RNEXT
                                                      0,   # PNEXT
                                                      0,   # TLEN
                                                      out_seq,   # SEQ
                                                      "*"]))))   # QUAL
            sam_writer.write("\n")

        seq_writer.close()
        fasta_writer.close()
        sam_writer.close()
"""
