import argparse
from os.path import isfile
import numpy as np
import pandas as pd
from multiprocessing import Pool
from logzero import logger
from BITS.scheduler import Scheduler
from BITS.utils import run_command, debug_mode, load_pickle, save_pickle
from .scheduler_args import add_scheduler_args

dir_name = "datruf"
tr_intervals_pkl = f"{dir_name}/tr_intervals.pkl"
alignments_pkl = f"{dir_name}/alignments.pkl"


def _find_units(read_id):
    # main subroutine


def _find_units_mult(read_ids):
    return [find_units(read_id) for read_id in read_ids]


def find_units(start_dbid, end_dbid, n_core):   # entry point of distributed computation   # TODO: how to pass the args without main and argparser?

    # Load dump data
    tr_intervals_all = (load_pickle(tr_intervals_pkl) \
                        .pipe(lambda df: start_dbid <= df["dbid"] and df["dbid"] <= end_dbid, axis=1))
    alignments_all = (load_pickle(alignments_pkl) \
                      .pipe(lambda df: start_dbid <= df["dbid"] and df["dbid"] <= end_dbid, axis=1))

    read_ids = sorted(set(tr_intervals_all["dbid"]))

    results = pd.DataFrame()
    units = pd.DataFrame()
    if n_core == 1:
        _find_units_mult(read_ids, tr_intervals_all, alignments_all)
    else:
        # split list
        n_per_core = -(-len(read_ids) // n_core)
        list_args = [(read_ids[i * n_per_core : (i + 1) * n_per_core],
                      tr_intervals_all,
                      alignments_all)
                     for i in range(n_core)]
        with Pool(n_core) as pool:
            for ret in pool.starmap(find_units_mult, list_args):
                # add results
                pass
    return (results, units)


def generate_dumps(args):
    # Extract data from DBdump's and LAdump's output
    save_pickle(pd.DataFrame([x.split('\t') for x in run_command(f"DBdump -r -h -mtan {args.db_file} | awk '$1 == \"R\" {{dbid = $2}} $1 == \"T0\" && {args.start_dbid} <= dbid && dbid <= {args.end_dbid} && $2 > 0 {{for (i = 1; i <= $2; i++) printf(\"%s\\t%s\\t%s\\n\", dbid, $(2 * i + 1), $(2 * i + 2))}}'").strip().split('\n')],
                             columns=("dbid", "start", "end")),
                tr_intervals_pkl)
    save_pickle(pd.DataFrame([x.split('\t') for x in run_command(f"LAdump -c {args.db_file} {args.las_file} | awk '$1 == \"P\" {{dbid = $2}} $1 == \"C\" && {args.start_dbid} <= dbid && dbid <= {args.end_dbid} {{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\", dbid, $2, $3, $4, $5)}}'").strip().split('\n')],
                             columns=("dbid", "abpos", "aepos", "bbpos", "bepos")),
                alignments_pkl)


def main():
    args = load_args()
    if args.end_id <= 0:   # Set <end_id> as the last read if not specified
        args.end_id = int(run_command(f"DBdump {args.db_file} | awk 'NR == 1 {{print $3}}'").strip())
    run_command(f"mkdir -p {dir_name}; rm -f {dir_name}/*")

    generate_dumps(args)

    if args.job_scheduler is None:
        find_units(args.start_dbid, args.end_dbid, args.n_core)
    else:
        n_per_job = -(-(args.end_dbid - args.start_dbid + 1) // n_distribute)
        n_digit = int(np.log10(args.n_distribute) + 1)
        s = Scheduler(args.job_scheduler,
                      args.submit_command,
                      args.queue_name,
                      f"{dir_name}/log")
        jids = []
        for i in range(args.n_distribute):
            index = str(i + 1).zfill(n_digit)
            jids.append(s.submit(' '.join([f"python -m dacembler.datruf",
                                           f"-s {args.start_dbid + i * n_per_job}",
                                           f"-e {min([args.start_dbid + (i + 1) * n_per_job - 1, args.end_dbid])}",
                                           f"-n {args.n_core}",
                                           f"{'-D' if args.debug_mode else ''}",
                                           f"--index {index}"
                                           f"{args.db_file}",
                                           f"{args.las_file}"]),
                                 f"{dir_name}/run_datruf.sh.{index}",
                                 job_name="datruf_distribute",
                                 n_core=args.n_core))

        # Merge the results
        logger.info("Waiting for all jobs to be finished...")
        s.submit('\n'.join([f"cat {dir_name}/datruf_results.* > {dir_name}/datruf_results.cat",
                            f"cat {dir_name}/datruf_units.* > {dir_name}/datruf_units.cat",
                            f"awk -F'\\t' 'BEGIN {{count = 0}} NR == 1 {{print $0}} $1 != \"\" {{printf count; for (i = 2; i <= NF; i++) {{printf \"\\t\" $i}}; print \"\"; count++}}' {dir_name}/datruf_results.cat > datruf_results",
                            f"awk -F'\\t' 'BEGIN {{count = 0}} NR == 1 {{print $0}} $1 != \"\" {{printf count; for (i = 2; i <= NF; i++) {{printf \"\\t\" $i}}; print \"\"; count++}}' {dir_name}/datruf_units.cat > datruf_units"]),
                 f"{dir_name}/finalize_datruf.sh",
                 job_name="datruf_finalize",
                 depend=jids,
                 wait=True)


def main_distribute():
    args = load_args()
    find_units(args.start_dbid, args.end_dbid, args.n_core)


def load_args():
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

    p.add_argument("--index",
                   type=int,
                   help="Used for distributed computation.")

    add_scheduler_args(p)
    args = p.parse_args()
    debug_mode(args.debug_mode)
    if args.n_distribute > 1:
        assert args.job_scheduler is not None, "job_scheduler must be specified"
        assert args.submit_command is not None, "submit_command must be specified"
    return args


if __name__ == "__main__":
    main_distribute()
