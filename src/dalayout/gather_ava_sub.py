import pandas as pd
from BITS.utils import run_command, load_pickle, save_pickle

def main():
    save_pickle(pd.concat([load_pickle(fname)
                           for fname in run_command(f"find . -name 'overlaps.*.pkl'").strip().split('\n')]) \
                .sort_values(by="strand") \
                .sort_values(by="read_j", kind="mergesort") \
                .sort_values(by="read_i", kind="mergesort") \
                .reset_index(drop=True),
                "overlaps.pkl")
