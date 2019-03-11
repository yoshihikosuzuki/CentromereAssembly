from os.path import isfile
import numpy as np
import pandas as pd
from multiprocessing import Pool
from logzero import logger
from BITS.utils import run_command, debug_mode


def _find_units(read_id):
    # main subroutine


def _find_units_mult(read_ids):
    return [find_units(read_id) for read_id in read_ids]


def find_units(tr_intervals_all, alignments_all, n_core):   # entry point of distributed computation   # TODO: how to pass the args without main and argparser?
    read_ids = set(tr_intervals_all["dbid"])

    if n_core == 1:
        _find_units_mult(read_ids, tr_intervals_all, alignments_all)
    else:
        # split list
        n_per_core =  -(-len(read_ids) // n_core)
        list_args = [(read_ids[i * n_per_core : (i + 1) * n_per_core],
                      tr_intervals_all,
                      alignments_all)
                     for i in range(n_core)]
        with Pool(n_core) as pool:
            for ret in pool.starmap(find_units_mult, list_args):
                # add results


def check_end_id(args):
    # Set <end_id> as the last read if not specified
    if args.end_id <= 0:
        args.end_id = int(run_command(f"DBdump {args.db_file} | awk 'NR == 1 {{print $3}}'").strip())


def load_dumps(args):
    # Extract data from DBdump's and LAdump's output
    # NOTE: start/end read IDs are both inclusive
    tr_intervals_all = pd.DataFrame([x.split('\t') for x in run_command(f"DBdump -r -h -mtan {args.db_file} | awk '$1 == \"R\" {{dbid = $2}} $1 == \"T0\" && {args.start_id} <= dbid && dbid <= {args.end_id} && $2 > 0 {{for (i = 1; i <= $2; i++) printf(\"%s\\t%s\\t%s\\n\", dbid, $(2 * i + 1), $(2 * i + 2))}}'").strip().split('\n')],
                                    columns=("dbid", "start", "end"))
    alignments_all = pd.DataFrame([x.split('\t') for x in run_command(f"LAdump -c {args.db_file} {args.las_file} | awk '$1 == \"P\" {{dbid = $2}} $1 == \"C\" && {args.start_id} <= dbid && dbid <= {args.end_id} {{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\", dbid, $2, $3, $4, $5)}}'").strip().split('\n')],
                                  columns=("dbid", "abpos", "aepos", "bbpos", "bepos"))
    return (tr_intervals_all, alignments_all)


def main():
    args = load_args()
    check_end_id(args)
    tr_intervals_all, alignments_all = load_dumps(args)
    find_units(tr_intervals_all, alignments_all, args.n_core)


def load_args():
    import argparse
    p = argparse.ArgumentParser(description="Find tandem repeat units from PacBio reads.")

    p.add_argument("db_file",
                   help="DAZZ_DB file.")

    p.add_argument("las_file",
                   help="Output of TANmask of the modified DAMASTER package.")

    p.add_argument("-s",
                   "--start_dbid",
                   type=int,
                   default=1,
                   help="Read ID of DAZZ_DB from which datruf's computation starts. [1]")

    p.add_argument("-e",
                   "--end_dbid",
                   type=int,
                   default=-1,
                   help="Read ID of DAZZ_DB at which datruf's computation ends. -1 means the last read. [-1]")

    p.add_argument("-n",
                   "--n_core",
                   type=int,
                   default=1,
                   help="Degree of parallelization. [1]")

    p.add_argument("-p",
                   "--n_distribute",
                   type=int,
                   default=1,
                   help="Degree of parallelization. [1]")

    p.add_argument("-D",
                   "--debug_mode",
                   action="store_true",
                   default=False,
                   help="Show debug messages. [False]")

    args = p.parse_args()
    debug_mode(args.debug_mode)
    return args


if __name__ == "__main__":
    main()
