import logging
import logzero
from logzero import logger
import pandas as pd
from BITS.utils import load_pickle, save_pickle
from .graph import Overlap


def main():
    args = load_args()

    # Calculate all-vs-all read alignments
    o = Overlap(pd.read_csv(args.reads_fname, sep='\t', index_col=0),
                load_pickle(args.encodings_fname),
                varvec_colname=args.varvec_colname)

    o.ava_read_alignment_distribute(args.n_distribute,
                                    args.n_core,
                                    args.job_scheduler,
                                    args.submit_command,
                                    args.queue_or_partition)


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=(""))

    parser.add_argument(
        "-r",
        "--reads_fname",
        type=str,
        default="tr_reads",
        help=("TR reads. [tr_reads]"))

    parser.add_argument(
        "-e",
        "--encodings_fname",
        type=str,
        default="encodings.pkl",
        help=("Encodings with variant vector units. [encodings.pkl]"))

    parser.add_argument(
        "-v",
        "--varvec_colname",
        type=str,
        default="var_vec_global0.15",
        help=("Column name to be used for alignemt. [var_vec_global0.15]"))

    parser.add_argument(
        "-p",
        "--n_distribute",
        type=int,
        default=1,
        help=("Degree of job distribution. [1]"))

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
