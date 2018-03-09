import argparse
import os
import re
import numpy as np
from sklearn.neighbors import KernelDensity
from interval import interval
import edlib

from datruf_utils import run_command


class Peaks:
    def __init__(self, peaks_dir, peak_fname_prefix):
        # NOTE: this class will have following variables:
        # During peak detection:
        #    self.units = [(header, seq), ...]
        #    self.peaks = [peak_len, ...]
        #    self.peak_intervals = [interval(peak_start, peak_end), ...]
        # During star-position alignment:
        #    self.n_peaks = <int>
        #    self.peak_units = {peak_id: [(header, seq), ...], ...}

        self.peaks_dir = peaks_dir
        self.peak_fname_prefix = peak_fname_prefix

    def load_fasta(self, fasta_fname):
        units = {}
        header = ""
        seq = ""
        with open(fasta_fname, 'r') as f:
            flag_first = True
            for line in f:
                line = line.strip()
                if line[0] == '>':
                    if flag_first:
                        flag_first = False
                    else:
                        units[header] = seq.upper()

                    header = line[1:]
                    if header in units:
                        print("[CAUTION] Header duplication!")
                    seq = ""
                else:
                    seq += line

            units[header] = seq.upper()

        return list(units.items())   # [(header, seq), ...]

    def load_units(self, fasta_fname):
        self.units = self.load_fasta(fasta_fname)

    def estimate_peaks(self, min_len=50, max_len=450, band_width=10):   # TODO: how to determine min/max_len? peak + density threshold?
        unit_lens = np.array([len(x[1]) for x in self.units])

        ul = unit_lens[(min_len < unit_lens) & (unit_lens < max_len)]
        ls = np.linspace(min_len, max_len, max_len - min_len + 1, dtype=int)

        kde = KernelDensity(kernel='gaussian', bandwidth=band_width).fit(ul.reshape(-1, 1))
        dens = np.exp(kde.score_samples(ls.reshape(-1, 1)))

        # Peak detection by naive sweep line
        self.peaks = []
        for i in range(1, len(dens) - 1):
            if (dens[i - 1] < dens[i]) and (dens[i] > dens[i + 1]):
                self.peaks.append(ls[i])
                print("[INFO] peak detected:", ls[i], dens[i])

    def calculate_peak_intervals(self, deviation=0.03):   # TODO: change deviation to percent of confidence interval
        # TODO: merge multiple neighboring peaks

        # Intervals around each peak unit length
        # Units will be collected from these intervals
        self.peak_intervals = interval()
        for x in self.peaks:
            self.peak_intervals |= interval[-(- x * (1. - deviation) // 1),
                                            int(x * (1. + deviation))]

    def output_peak_units(self):
        run_command("mkdir -p %s" % self.peaks_dir)

        writers = [open("%s/%s.%d.fasta"
                        % (self.peaks_dir,
                           self.peak_fname_prefix,
                           i),
                        'w')
                   for i in range(1, len(self.peaks) + 1)]

        for header, seq in self.units:
            for i, intvl in enumerate(self.peak_intervals.components):
                if len(seq) in intvl:
                    writers[i].write(">%s\n%s\n" % (header, seq))
                    break

        for i in range(len(writers)):
            writers[i].close()

    def load_peaks(self):
        command = "ls -l %s/%s.[0-9*].fasta | wc -l" % (self.peaks_dir,
                                                        self.peak_fname_prefix)
        self.n_peaks = int(run_command(command))

        self.peak_units = {}
        for peak_id in range(1, self.n_peaks + 1):   # peak ID must be 1, 2, 3, ...
            self.peak_units[peak_id] = self.load_fasta("%s/%s.%d.fasta"
                                                       % (self.peaks_dir,
                                                          self.peak_fname_prefix,
                                                          peak_id))

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


def main():
    args = load_args()

    # Peaks are yet detected
    if args.unit_fasta_fname is not None:
        peaks = Peaks(args.peaks_dir, args.peak_fname_prefix)
        peaks.load_units(args.unit_fasta_fname)
        peaks.estimate_peaks()
        peaks.calculate_peak_intervals()
        peaks.output_peak_units()
        del peaks

    if not os.path.isdir(args.peaks_dir):
        print("[ERROR] peak directory does not exist.")
        exit(1)

    peaks = Peaks(args.peaks_dir, args.peak_fname_prefix)
    peaks.load_peaks()
    peaks.align_start_positions()


def load_args():
    parser = argparse.ArgumentParser(
        description=("Run dacenter for peak units detection."))

    parser.add_argument(
        "--unit_fasta_fname",
        default=None,
        help=("Input fasta file of the unit sequences [None]"))

    parser.add_argument(
        "--peaks_dir",
        default="peaks",
        help=("Given --unit_fasta_fname, dacenter will find peaks from units "
              "and output units around the peaks into this directory. If not, "
              "dacenter will load peak units from this directory and perform "
              "start-position alignment of the units for each peak. [peaks]"))

    parser.add_argument(
        "--peak_fname_prefix",
        default="peak",
        help=("Prefix for output/input files of the peak units [peak]"))

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()
