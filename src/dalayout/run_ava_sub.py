import os.path
import logging
import logzero
from logzero import logger
import pandas as pd
from BITS.utils import load_pickle, save_pickle
from .encode import encode_reads, detect_variants


def main():
    args = load_args()

    overlap_obj = load_pickle(args.ovlp_pkl)
    list_pairs = load_pickle(args.list_pairs_pkl)
    overlap_obj._ava_read_alignment(list_pairs, args.n_core)
    save_pickle(overlap_obj.overlaps, args.out_pkl)


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=(""))

    parser.add_argument(
        "ovlp_pkl",
        type=str,
        help=("overlap_obj.pkl"))

    parser.add_argument(
        "list_pairs_pkl",
        type=str,
        help=("list_pairs_XX.pkl"))

    parser.add_argument(
        "out_pkl",
        type=str,
        help=("overlap"))

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
