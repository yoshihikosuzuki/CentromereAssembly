import logging
import logzero
from logzero import logger
import numpy as np
import pandas as pd
from BITS.utils import load_pickle, save_pickle, run_command


def main():
    args = load_args()

    c = load_pickle(args.clustering_obj_pkl)
    c._calc_dist_mat(load_pickle(args.rows_pkl), args.n_core)
    save_pickle(c._calc_dist_mat(load_pickle(args.rows_pkl), args.n_core), args.out_pkl)


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=("Perform layout of reads based on the representative units given."))

    parser.add_argument(
        "clustering_obj_pkl",
        type=str,
        help=("Clustering pickle object file."))

    parser.add_argument(
        "rows_pkl",
        type=str,
        help=("Row indices pickle file."))

    parser.add_argument(
        "out_pkl",
        type=str,
        help=("Output list of dist arrays pickle file."))

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
        pd.set_option('expand_frame_repr', False)   # show entire dataframe
    else:
        logzero.loglevel(logging.INFO)
    del args.debug_mode

    return args


if __name__ == "__main__":
    main()
