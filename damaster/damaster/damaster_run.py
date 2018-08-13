import argparse
import copy
import pickle
import logging
import logzero
from logzero import logger

from .damaster_core import PeaksFinder
from .damaster_io import load_peaks


def main():
    args = load_args()

    # Detect peaks in the unit length distribution
    finder = PeaksFinder(args.unit_fasta)
    finder.run()

    # Keep only Peak instances and discard the others
    peaks = copy.deepcopy(finder.peaks)
    del finder

    # Define a set of representative monomers from somewhat homogeneous TRs
    for peak in peaks:
        peak.construct_master_units(min_n_units=args.min_n_units,
                                    n_core=args.n_core)

    # NOTE: if you want to update data in peaks, use the pickle as much as possible

    # Output the peaks as pickle
    # Now output the list itself instead of each peak for simplicity
    with open("peaks.pkl", 'wb') as f:
        pickle.dump(peaks, f)


def load_args():
    parser = argparse.ArgumentParser(
        description=("Construct master units and representative units."))

    parser.add_argument(
        "unit_fasta",
        help=("Input fasta file of the unit sequences reported by datruf. "
              "[datruf_units.fasta]"))

    parser.add_argument(
        "-m",
        "--min_n_units",
        type=int,
        default=10,
        help=("Minimum number of units in a single TR whose consensus will be"
              "taken. [10]"))

    parser.add_argument(
        "-n",
        "--n_core",
        type=int,
        default=1,
        help=("Degree of parallelization. [1]"))

    parser.add_argument(
        "-D",
        "--debug_mode",
        action="store_true",
        default=False,
        help=("Run in debug mode. [False]"))

    args = parser.parse_args()
    if args.debug_mode:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    del args.debug_mode

    return args


if __name__ == "__main__":
    main()
