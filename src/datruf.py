import os.parh.isfile
import numpy as np
import pandas as pd
from multiprocessing import Pool
from logzero import logger
from BITS.utils import run_command, debug_mode


def _find_units(read_id):
    # main subroutine


def _find_units_mult(read_ids):
    return [find_units(read_id) for read_id in read_ids]


def find_units(args):   # entry point of distributed computation   # TODO: how to pass the args without main and argparser?
    # load data

    # split list
    all_read_ids = list(range(args.start_id, args.end_id + 1))
    n_per_core =  -(-len(all_read_ids) // args.n_core)
    list_read_ids = [all_read_ids[i * n_per_core : (i + 1) * n_per_core]   # TODO: add any arguments needed
                     for i in range(args.n_core)]
    
    with Pool(n_core) as pool:
        for ret in pool.starmap(find_units_mult, list_read_ids):
            # add results


def check_input(args, dbdump="datander_dbdump", ladump="datander_ladump"):
    # TODO: convert dump files into DF-friendly format
    if not os.path.isfile(dbdump):
        logger.info("DBdump file does not exist. Generating.")
        run_command(f"DBdump -r -h -mtan {args.db_file} | awk '$1 == \"R\" {{dbid = $2}} $1 == \"T0\" && $2 > 0 {{for (i = 0; i < )}}' > {dbdump}")
    if not os.path.isfile(ladump):
        logger.info("LAdump file does not exist. Generating.")
        run_command(f"LAdump -c {args.db_file} {args.las_file} > {ladump}")


def main():
    args = load_args()
    check_input(args)
    if args.end_id <= 0:   # last read
        args.end_id = int(run_command(f"DBdump {args.db_file} | awk 'NR == 1 {{print $3}}'").strip())
    find_units(args)


def load_args():
    import argparse
    p = argparse.ArgumentParser(description="Find tandem repeat units from PacBio reads.")

    p.add_argument("db_file",
                   help="DAZZ_DB file.")

    p.add_argument("las_file",
                   help="Output of TANmask of the modified DAMASTER package.")

    p.add_argument("-d",
                   "--dbdump",
                   type=str,
                   default="datander_dbdump",
                   help="DBdump file which will be automatically generated. [datander_dbdump]")

    p.add_argument("-l",
                   "--ladump",
                   type=str,
                   default="datander_ladump",
                   help="LAdump file which will be automatically generated. [datander_ladump]")

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

    p.add_argument("-D",
                   "--debug_mode",
                   action="store_true",
                   default=False,
                   help="Show debug messages. [False]")

    args = parser.parse_args()
    debug_mode(args.debug_mode)
    return args


if __name__ == "__main__":
    main()
