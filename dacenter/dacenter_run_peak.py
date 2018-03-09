import argparse
import os

from dacenter_peak import Peaks


def main():
    args = load_args()

    # Peaks are yet detected
    if args.unit_fasta_fname is not None:
        #if os.path.isdir(args.peaks_dir):
        #    print("[ERROR] peak directory already exists. You must remove it"
        #          "before running this.")
        #    exit(1)

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
