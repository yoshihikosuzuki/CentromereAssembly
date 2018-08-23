import os
import logging
import logzero
from logzero import logger
import pandas as pd
from .peak import save_peaks, load_peaks
from .clustering import ClusteringSeqs


def main():
    args = load_args()
    pkl_fname = "peaks.pkl"

    if not args.from_scratch and os.path.isfile(pkl_fname):
        logger.info(f"Loading data from existing {pkl_fname}")
        peaks = load_peaks(pkl_fname)
    else:
        logger.info(f"Data will be stored to {pkl_fname} and re-used next time")

        # Detect peaks in the unit length distribution
        logger.info("Starting peak detection")
        from .peak import PeaksFinder
        peaks = PeaksFinder(args.units_fname).detect_peaks()
        save_peaks(peaks, pkl_fname)

    for i, peak in enumerate(peaks):
        if i != len(peaks) - 1:   # only the last peak   # NOTE: for debug
            continue

        # Generate intra-TR consensus units
        if args.from_scratch or not hasattr(peak, "cons_units"):
            logger.info("Starting taking intra-TR consensus")
            peak.take_intra_consensus(args.min_n_units, args.n_core)
            logger.info("Finished intra-TR consensus")
            save_peaks(peaks, pkl_fname)

        # Cluster the intra-TR consensus units
        if args.from_scratch or not hasattr(peak, "cl_master"):   # TODO: incorporate from_scratch and clustering
            peak.cl_master = ClusteringSeqs(peak.cons_units["sequence"])
        logger.info("Starting hierarchical clustering")
        peak.cl_master.cluster_hierarchical(n_core=args.n_core)
        logger.info("Finished hierarchical clustering")
        save_peaks(peaks, pkl_fname)

        # Generate master units
        peak.master_original = peak.cl_master.generate_consensus()
        logger.debug(f"\n{peak.master_original}")
        peak.filter_master_units()
        logger.info("Finished cluster consensus")
        save_peaks(peaks, pkl_fname)

        #peak.construct_repr_units(n_core=args.n_core)


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=("Construct master units and representative units."))

    parser.add_argument(
        "-u",
        "--units_fname",
        type=str,
        default="datruf_units",
        help=("Input file of the unit sequences reported by datruf. [datruf_units]"))

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
        pd.set_option('expand_frame_repr', False)   # to show full df
    else:
        logzero.loglevel(logging.INFO)
    del args.debug_mode

    return args


if __name__ == "__main__":
    main()
