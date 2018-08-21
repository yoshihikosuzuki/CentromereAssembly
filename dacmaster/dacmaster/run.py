import argparse
import copy
import pickle
import logging
import logzero
from logzero import logger
import pandas as pd

from .peak import Peak, PeaksFinder
from .clustering import ClusteringSeqs


def load_precomputed(precomputed):
    # <precomputed["peaks"]> corresponds to <PeaksFinder.peaks>
    # Re-create Peak instances so that any modifications on Peak class would not affect loading
    # If any of these members are modified, re-claculate from the begining
    peaks = [Peak(p.peak_id, p.N, p.unit_len, p.density, p.min_len, p.max_len, p.raw_units)
             for p in precomputed["peaks"]]

    # These data below take time to calculate, thus re-use if previous one is available
    if "units_consensus" in precomputed:
        for i, data in enumerate(precomputed["units_consensus"]):
            if data is not None:
                setattr(peaks[i], "units_consensus", data)
                setattr(peaks[i], "clustering", ClusteringSeqs(peaks[i].units_consensus["sequence"]))
    if "dist_matrix" in precomputed:   # TODO: rename?
        for i, data in enumerate(precomputed["dist_matrix"]):
            if data is not None:
                setattr(peaks[i].clustering, "dist_matrix", data)
    if "hc_result_precomputed" in precomputed:
        for i, data in enumerate(precomputed["hc_result_precomputed"]):
            if data is not None:
                setattr(peaks[i].clustering, "hc_result_precomputed", data)

    return peaks


def main():
    args = load_args()

    # You can skip redundant heavy calculation by specifying precomputed pickle
    if args.precomputed_pkl is not None:
        logger.info(f"Loading precomputed data from {args.precomputed_pkl}")
        with open(args.precomputed_pkl, 'rb') as f:
            precomputed = pickle.load(f)
    else:
        precomputed = {}
        setattr(args, "precomputed_pkl", "precomputed.pkl")
        logger.info(f"As well as peaks.pkl, some computationally heavy data are "
                    f"stored to {args.precomputed_pkl} for the next time.")

    ## -------------------------------------- ##
    ## Step 1. Detection of peak unit lengths ##
    ## -------------------------------------- ##

    if "peaks" in precomputed:
        peaks = load_precomputed(precomputed)
    else:
        logger.info("No precomputed Peak instances. Staring peak detection")
        # Detect peaks in the unit length distribution
        finder = PeaksFinder(args.units_fname)
        finder.run()

        # Keep only Peak instances and discard the others
        peaks = copy.deepcopy(finder.peaks)
        del finder

        # Update precomputed data
        # Regardless of the following computation, save it here
        precomputed["peaks"] = peaks
        with open(args.precomputed_pkl, 'wb') as f:
            pickle.dump(precomputed, f)
        logger.info(f"Saved Peak instances to {args.precomputed_pkl}")

    for i, peak in enumerate(peaks):
        if i != len(peaks) - 1:   # only the last peak   # NOTE: for debug
            continue

        ## ------------------------------------ ##
        ## Step 2. Construction of master units ##
        ## ------------------------------------ ##

        # Cluster intra-TR consensus unit sequences
        peak.construct_master_units(args.min_n_units, args.n_core)

        ## -------------------------------------------- ##
        ## Step 3. Construction of representative units ##
        ## -------------------------------------------- ##

        #peak.construct_repr_units(n_core=args.n_core)

    # Update precomputed data (even when there is actually no update!)
    precomputed["units_consensus"] = [peak.units_consensus
                                      if hasattr(peak, "units_consensus")
                                      else None
                                      for peak in peaks]
    precomputed["dist_matrix"] = [peak.clustering.dist_matrix
                                  if hasattr(peak, "clustering")
                                  and hasattr(peak.clustering, "dist_matrix")
                                  else None
                                  for peak in peaks]
    precomputed["hc_result_precomputed"] = [peak.clustering.hc_result_precomputed
                                            if hasattr(peak, "clustering")
                                            and hasattr(peak.clustering, "hc_result_precomputed")
                                            and len(peak.clustering.hc_result_precomputed) > 0
                                            else None
                                            for peak in peaks]
    with open(args.precomputed_pkl, 'wb') as f:
        pickle.dump(precomputed, f)
    logger.info(f"Saved heavy data in Peak instances to {args.precomputed_pkl}")

    # Output the whole result
    with open("peaks.pkl", 'wb') as f:
        pickle.dump(peaks, f)


def load_args():
    parser = argparse.ArgumentParser(
        description=("Construct master units and representative units."))

    parser.add_argument(
        "units_fname",
        help=("Input file of the unit sequences reported by datruf."))

    parser.add_argument(
        "-m",
        "--min_n_units",
        type=int,
        default=10,
        help=("Minimum number of units in a single TR whose consensus will be"
              "taken. [10]"))

    parser.add_argument(
        "-p",
        "--precomputed_pkl",
        type=str,
        default=None,
        help=("Pickle file keeping precomputed variables. [None]"))

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
        pd.set_option('expand_frame_repr', False)   # to show full df
    else:
        logzero.loglevel(logging.INFO)
    del args.debug_mode

    return args


if __name__ == "__main__":
    main()
