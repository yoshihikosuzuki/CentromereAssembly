import os
import argparse
import numpy as np
from logzero import logger

from BITS import utils


def check_dump(args):
    """
    Generate dump files beforehand instead of in each distributed job.
    """

    if not os.path.isfile(args.dbdump):
        command = f"DBdump -r -h -mtan {args.db_file} > {args.dbdump}"
        utils.run_command(command)

    if not os.path.isfile(args.ladump):
        command = f"LAdump -c {args.db_file} {args.las_file} > {args.ladump}"
        utils.run_command(command)


def main():
    args = load_args()

    check_dump(args)

    command = f"DBdump {args.db_file} | awk 'NR == 1 {{print $3}}'"
    max_dbid = int(utils.run_command(command).strip())

    if args.start_dbid < 1:
        args.start_dbid = 1
    if args.end_dbid < 1 or args.end_dbid > max_dbid:
        args.end_dbid = max_dbid

    n_dbid_part = -(-args.end_dbid // args.n_distribute)
    n_digit = int(np.log10(args.n_distribute) + 1)

    for i in range(args.n_distribute):
        idx = str(i + 1).zfill(n_digit)

        script_fname = f"run_datruf.{args.job_scheduler}.{idx}"
        with open(script_fname, 'w') as f:
            start = args.start_dbid + i * n_dbid_part
            end = args.start_dbid + (i + 1) * n_dbid_part - 1
            script = (f"datruf_run.py -s {start} -e {end} -n {args.n_core} "
                      f"-m {args.out_main_fname}.{idx} -u {args.out_units_fname}.{idx} "
                      f"{'--only_interval' if args.only_interval else ''} "
                      f"{args.db_file} {args.las_file}")
            if args.job_scheduler == "sge":
                script = utils.sge_nize(script,
                                        job_name="run_datruf",
                                        n_core=args.n_core,
                                        sync=False)
            elif args.job_scheduler == "slurm":
                script = utils.slurm_nize(script,
                                          job_name="run_datruf",
                                          n_core=args.n_core,
                                          mem_per_cpu=40000,
                                          wait=False)
            f.write(script)

        command = f"{args.submit_job} {script_fname}"
        utils.run_command(command)

    # Generate a script that should be ran after the datruf calculation finishes
    with open("finalize_datruf.sh", 'w') as f:
        f.write(
f"""cat {args.out_units_fname}.* > {args.out_units_fname}
cat {args.out_main_fname}.* > {args.out_main_fname}.cat
awk -F'\\t' 'NR == 1 {{print $0}} $1 != \"\" {{print $0}}' {args.out_main_fname}.cat > {args.out_main_fname}
rm {args.out_units_fname}.*; rm {args.out_main_fname}.*
""")

    logger.info("Run `$ bash finalize_datruf.sh` after finishing all jobs")


def load_args():
    parser = argparse.ArgumentParser(description="Run datruf with many reads using job scheduler")

    parser.add_argument(
        "db_file",
        help="DAZZ_DB file")

    parser.add_argument(
        "las_file",
        help="output of TANmask")

    parser.add_argument(
        "-d",
        "--dbdump",
        type=str,
        default="datander_dbdump",
        help=("Output of `DBdump -r -h -mtan <db_file>`. This will be "
              "automatically generated if not exist. [datander_dbdump]"))

    parser.add_argument(
        "-l",
        "--ladump",
        type=str,
        default="datander_ladump",
        help=("Output of `LAdump -c <db_file> <las_file>`. This will be "
              "automatically generated if not exist. [datander_ladump]"))

    parser.add_argument(
        "-s",
        "--start_dbid",
        type=int,
        default=1,
        help=("Start read ID, which is used in DAZZ_DB. Set <= 1 to start "
              "from the first read. [1]"))

    parser.add_argument(
        "-e",
        "--end_dbid",
        type=int,
        default=-1,
        help=("End read ID. Set < 1 to end at the last read. [-1]"))

    parser.add_argument(
        "-n",
        "--n_core",
        type=int,
        default=1,
        help=("Degree of parallelization in each distributed job. [1]"))

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
        "--submit_job",
        type=str,
        default="qsub",
        help="Command name to submit a job with the specified scheduler. [qsub]")

    parser.add_argument(
        "-m",
        "--out_main_fname",
        type=str,
        default="datruf_result",
        help=("Write main results to this file. [datruf_result]"))

    parser.add_argument(
        "-u",
        "--out_units_fname",
        type=str,
        default="datruf_units.fasta",
        help=("Write unit sequences to this file. [datruf_units.fasta]"))

    parser.add_argument(
        "--only_interval",
        action="store_true",
        default=False,
        help=("Stop calculation just after obtaining TR intervals. Also "
              "filtering by CV of the unit lengths is not applied. [False]"))

    return parser.parse_args()


if __name__ == "__main__":
    main()
