import os
import sys
import pickle
import numpy as np
import pandas as pd
from logzero import logger


def load_unit_fasta(fasta_fname):
    """
    Return unit sequences in a fasta file as a pandas dataframe with the following columns:

    <index> <read header> <read ID> <path ID> <unit ID> <start pos> <end pos> <length> <sequence>
    """

    units = {}
    header = ""
    index = 0
    with open(fasta_fname, 'r') as f:
        flag_first = True
        for line in f:
            line = line.strip()
            if line[0] == '>':
                if flag_first:
                    flag_first = False
                else:
                    units[index] = (header,
                                    read_id,
                                    path_id,
                                    unit_id,
                                    start,
                                    end,
                                    read_len,
                                    seq.upper())
                    index += 1

                header = line[1:]
                first, second, third = header.split('/')
                read_id, path_id = list(map(int, first.split('-')))
                unit_id = int(second)
                start, end = list(map(int, third.split(' ')[0].split('_')))
                read_len = end - start
                seq = ""
            else:
                seq += line

        units[index] = (header,
                        read_id,
                        path_id,
                        unit_id,
                        start,
                        end,
                        read_len,
                        seq.upper())

    return pd.DataFrame.from_dict(units,
                                  orient="index",
                                  columns=("header",
                                           "read_id",
                                           "path_id",
                                           "unit_id",
                                           "start",
                                           "end",
                                           "length",
                                           "sequence"))


def load_peaks(peaks_dir):
    """
    Load a list of Peak instances.
    """

    pkl_fname = os.path.join(peaks_dir, "peaks.pkl")
    if not os.path.isfile(pkl_fname):
        logger.error(f"{pkl_fname} does not exist. Abort.")
        sys.exit(1)

    with open(pkl_fname, 'rb') as f:
        peaks = pickle.load(f)

    if len(peaks) == 0:
        logger.error("No peak pikles were detected.")
    else:
        logger.info(f"{len(peaks)} peaks were loaded.")

    return peaks
