import pandas as pd
from BITS.utils import run_command, load_pickle, save_pickle


def main():
    args = load_args()
    save_pickle(pd.concat([load_pickle(fname)
                           for fname in run_command(f"find {args.dir_name} -name 'overlaps.*.pkl'").strip().split('\n')]) \
                .sort_values(by="strand") \
                .sort_values(by="read_j", kind="mergesort") \
                .sort_values(by="read_i", kind="mergesort") \
                .reset_index(drop=True),
                args.out_pkl_fname)


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=("Merge results of ava read overlaps."))

    parser.add_argument(
        "-d",
        "--dir_name",
        type=str,
        default=".",
        help=("Directory name in which overlaps.*.pkl exist. [.]"))

    parser.add_argument(
        "-o",
        "--out_pkl_fname",
        type=str,
        default="overlaps.pkl",
        help=("Output pickle file for overlaps. [overlaps.pkl]"))

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
