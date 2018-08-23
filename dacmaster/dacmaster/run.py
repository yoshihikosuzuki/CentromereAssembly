import os
import pickle
import logging
import logzero
from logzero import logger
import numpy as np
import pandas as pd
from .peak import Peak, PeaksFinder
from .clustering import ClusteringSeqs


# TODO: save whole Peak instances into a pickle, and load only necessary variables and create new instances

"""
@dataclass(repr=False, eq=False)
class Precomputed:
    n_peak    : InitVar[int]
    peaks     : list = field(init=False)
    cons_units: list = field(init=False)
    dist_mat  : list = field(init=False)
    hc_res_pre: list = field(init=False)

    def __post_init__(self):
        self.peaks = 
"""


def load_precomputed(precomputed):
    # <precomputed["peaks"]> corresponds to <PeaksFinder.peaks>
    # Re-create Peak instances so that any modifications on Peak class would not affect loading
    # If any of these members are modified, re-claculate from the begining
    peaks = [Peak(p.peak_id, p.N, p.unit_len, p.density, p.min_len, p.max_len, p.raw_units)
             for p in precomputed["peaks"]]

    # These data below take time to calculate, thus re-use if previous one is available
    if "cons_units" in precomputed:
        for i, data in enumerate(precomputed["cons_units"]):
            if data is not None:
                setattr(peaks[i], "cons_units", data)
                setattr(peaks[i], "cl_master", ClusteringSeqs(peaks[i].cons_units["sequence"]))
    if "dist_matrix" in precomputed:   # TODO: rename?
        for i, data in enumerate(precomputed["dist_matrix"]):
            if data is not None:
                setattr(peaks[i].cl_master, "dist_matrix", data)
    if "hc_result_precomputed" in precomputed:
        for i, data in enumerate(precomputed["hc_result_precomputed"]):
            if data is not None:
                setattr(peaks[i].cl_master, "hc_result_precomputed", data)

    return peaks


def main():
    args = load_args()

    # You can skip redundant heavy calculation by specifying precomputed pickle
    if os.path.isfile(args.precomputed_pkl):
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
        peaks = finder.peaks   # TODO: maybe should rollback to previous one

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
        #peak.construct_master_units(args.min_n_units, args.n_core)

        # Calculate intra-TR consensus unit for each TR
        if not hasattr(peak, "cons_units"):
            logger.info("Starting taking intra-TR consensus")
            peak.take_intra_consensus(args.min_n_units, args.n_core)
            logger.info("Finished intra-TR consensus")

        # Save data for each loop
        precomputed["cons_units"] = [peak.cons_units
                                     if hasattr(peak, "cons_units")
                                     else None
                                     for peak in peaks]
        with open(args.precomputed_pkl, 'wb') as f:
            pickle.dump(precomputed, f)
        logger.info(f"Saved consensus units to {args.precomputed_pkl}")

        # Cluster the intra-TR consensus units
        if not hasattr(peak, "cl_master"):
            peak.cl_master = ClusteringSeqs(peak.cons_units["sequence"])
        logger.info("Starting hierarchical clustering")
        peak.cl_master.cluster_hierarchical(n_core=args.n_core)
        logger.info("Finished hierarchical clustering")

        # Save data for each loop
        precomputed["dist_matrix"] = [peak.cl_master.dist_matrix
                                      if hasattr(peak, "cl_master")
                                      and hasattr(peak.cl_master, "dist_matrix")
                                      else None
                                      for peak in peaks]
        precomputed["hc_result_precomputed"] = [peak.cl_master.hc_result_precomputed
                                                if hasattr(peak, "cl_master")
                                                and hasattr(peak.cl_master, "hc_result_precomputed")
                                                and len(peak.cl_master.hc_result_precomputed) > 0
                                                else None
                                                for peak in peaks]
        with open(args.precomputed_pkl, 'wb') as f:
            pickle.dump(precomputed, f)
        logger.info(f"Saved clustering data to {args.precomputed_pkl}")

        # Generate master units
        peak.master_original = peak.cl_master.generate_consensus()
        logger.debug(f"\n{peak.master_original}")
        # remove strange sequences and redundancy
        # after that, start-position adjustment and strand adjustment of close seqs
        peak.filter_master_units()
        logger.info("Finished cluster consensus")

        ## -------------------------------------------- ##
        ## Step 3. Construction of representative units ##
        ## -------------------------------------------- ##

        #peak.construct_repr_units(n_core=args.n_core)

    with open(args.precomputed_pkl, 'wb') as f:
        pickle.dump(precomputed, f)
    logger.info(f"Saved heavy data in Peak instances to {args.precomputed_pkl}")

    # Output the whole result
    with open("peaks.pkl", 'wb') as f:
        pickle.dump(peaks, f)


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
        "-p",
        "--precomputed_pkl",
        type=str,
        default="precomputed.pkl",
        help=("Pickle file keeping precomputed variables. [precomputed.pkl]"))

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
