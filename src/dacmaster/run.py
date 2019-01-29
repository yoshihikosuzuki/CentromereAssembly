import os.path
import logging
import logzero
from logzero import logger
import pandas as pd
from .peak import Peaks, load_peaks
from BITS.utils import run_command


def main():
    args = load_args()

    # Generate reads DataFrame if needed
    if args.from_scratch or not os.path.isfile("reads"):
        logger.info(f"Generating reads DataFrame")
        run_command(f"DBshow -w10000000 {args.db_file} | "
                    f"awk -F'>' 'BEGIN {{print \"\tdbid\theader\tlength\tsequence\"; count = 1}} "
                    f"count % 2 == 1 {{header = $2}} "
                    f"count % 2 == 0 {{print (count / 2 - 1) \"\t\" (count / 2) \"\t\" "
                    f"header \"\t\" length($1) \"\t\" $1}} "
                    f"{{count++}}' > reads")
    assert os.path.getsize("reads") != 0, "The file 'reads' is empty!"

    # Detect peak unit lenghs
    if not args.from_scratch and os.path.isfile(args.out_pkl_fname):
        logger.info(f"Loading data from existing {args.out_pkl_fname}")
        peaks = load_peaks(args.out_pkl_fname)
    else:
        logger.info(f"Data will be stored to {args.out_pkl_fname} and re-used next time")

        # Detect peaks in the unit length distribution
        peaks = Peaks("reads", args.units_fname)
        peaks.detect_peaks()
        peaks.save(args.out_pkl_fname)

    # For each peak, calculate representative units
    for i, peak in enumerate(peaks.peaks):
        if i != len(peaks.peaks) - 1:   # only the last peak   # NOTE: for debug
            continue

        # Generate intra-TR consensus units
        if not hasattr(peak, "cons_units"):
            peak.take_intra_consensus(args.min_n_units, args.n_core)
            peaks.save(args.out_pkl_fname)

        # Cluster the intra-TR consensus units
        peak.cluster_cons_units(args.n_core)
        peaks.save(args.out_pkl_fname)

        # Generate master units
        peak.generate_master_units()
        peaks.save(args.out_pkl_fname)

        #peak.construct_repr_units(n_core=args.n_core)


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=("Construct master units and representative units."))

    parser.add_argument(
        "db_file",
        help=("DAZZ_DB file. Used for generating reads fasta and dbid_header."))

    parser.add_argument(
        "-u",
        "--units_fname",
        type=str,
        default="datruf_units",
        help=("Input file of the unit sequences reported by datruf. [datruf_units]"))

    parser.add_argument(
        "-o",
        "--out_pkl_fname",
        type=str,
        default="peaks.pkl",
        help=("Output pickle file. [peaks.pkl]"))

    parser.add_argument(
        "-m",
        "--min_n_units",
        type=int,
        default=10,
        help=("Minimum number of units in a single TR whose consensus will be taken. [10]"))

    parser.add_argument(
        "-n",
        "--n_core",
        type=int,
        default=1,
        help=("Degree of parallelization. [1]"))

    parser.add_argument(
        "-F",
        "--from_scratch",
        action="store_true",
        default=False,
        help=("Compute from scratch without loading existing peaks.pkl file. [False]"))

    parser.add_argument(
        "-D",
        "--debug_mode",
        action="store_true",
        default=False,
        help=("Run in debug mode. [False]"))

    args = parser.parse_args()
    if args.debug_mode:
        logzero.loglevel(logging.DEBUG)
        pd.set_option('expand_frame_repr', False)   # show entire dataframe
    else:
        logzero.loglevel(logging.INFO)
    del args.debug_mode

    return args


if __name__ == "__main__":
    main()
