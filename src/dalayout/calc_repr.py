import logging
import logzero
from logzero import logger
import numpy as np
import pandas as pd
from BITS.utils import load_pickle, save_pickle, run_command
from .encode import encode_reads, cut_unit_from_read
from dacmaster.clustering import ClusteringSeqs


def main():
    args = load_args()

    reads = pd.read_csv(args.tr_reads, sep='\t', index_col=0)
    encodings = load_pickle(args.encodings_fname)
    cover_rate = pd.read_csv(args.cover_rate, sep='\t', index_col=0)

    # Collect all synchronized raw units belonging to the peak_id and repr_id
    units = (encodings[encodings.apply(lambda df: (df["peak_id"] == args.peak_id
                                                   and df["repr_id"] == args.repr_id
                                                   and df["type"] != "boundary"
                                                   and cover_rate.loc[df["read_id"]]["cover_rate"] >= 0.5),
                                       axis=1)]
             .apply(lambda d: cut_unit_from_read(reads, d, False), axis=1))

    # Clustering
    c = ClusteringSeqs(units.reset_index(drop=True),
                       np.array(units.index),
                       cyclic=False,
                       rc=False)
    c.calc_dist_mat(args.n_core, args.n_distribute)
    c.cluster_hierarchical()
    c.generate_consensus()   # representative units
    save_pickle(c, "repr/clustering.{peak_id}.{repr_id}.pkl")


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=("Perform layout of reads based on the representative units given."))

    parser.add_argument(
        "peak_id",
        type=int,
        help=("Peak id"))

    parser.add_argument(
        "repr_id",
        type=int,
        help=("Master id"))

    parser.add_argument(
        "-e",
        "--encodings_fname",
        type=str,
        default="encodings.pkl",
        help=("Encodings by master units. [encodings.pkl]"))

    parser.add_argument(
        "-c",
        "--cover_rate",
        type=str,
        default="cover_rate",
        help=("Cover rate file by masters. [cover_rate]"))

    parser.add_argument(
        "-r",
        "--tr_reads",
        type=str,
        default="tr_reads",
        help=("TR reads. [tr_reads]"))

    parser.add_argument(
        "-n",
        "--n_core",
        type=int,
        default=1,
        help=("Degree of parallelization. [1]"))

    parser.add_argument(
        "-p",
        "--n_distribute",
        type=int,
        default=1,
        help=("Degree of parallelization in each distributed job. [1]"))

    parser.add_argument(
        "-j",
        "--job_scheduler",
        type=str,
        default="sge",
        help="Job scheduler name. ('sge' or 'slurm)' [sge]")

    parser.add_argument(
        "-c",
        "--submit_command",
        type=str,
        default="qsub",
        help="Command name to submit a job with the specified scheduler. [qsub]")

    parser.add_argument(
        "-q",
        "--queue_or_partition",
        type=str,
        default=None,
        help="Name of queue (SGE) or partition (SLURM) to which jobs are submitted. [None]")

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
