import argparse
import numpy as np
from logzero import logger

from .run import Runner
from BITS.utils import run_command, sge_nize, slurm_nize


def main():
    args = load_args()

    # Use Runner instance just for some preparation of data
    r = Runner(args)   # args are now variables of <r>
    r._check_dump()

    n_dbid_part = -(-r.end_dbid // args.n_distribute)
    n_digit = int(np.log10(args.n_distribute) + 1)

    # Submit scripts each of which runs datruf_run.py
    for i in range(args.n_distribute):
        index = str(i + 1).zfill(n_digit)
        script_fname = f"run_datruf.{args.job_scheduler}.{index}"

        with open(script_fname, 'w') as f:
            script = ' '.join([f"datruf_run.py",
                               f"-s {r.start_dbid + i * n_dbid_part}",
                               f"-e {r.start_dbid + (i + 1) * n_dbid_part - 1}",
                               f"-m {args.out_main_fname}.{index}",
                               f"-u {args.out_units_fname}.{index}",
                               f"{'--only_interval' if args.only_interval else ''}",
                               f"{'-D' if args.debug_mode else ''}",
                               f"-n {args.n_core}",
                               f"{args.db_file}",
                               f"{args.las_file}"])

            if args.job_scheduler == "sge":
                script = sge_nize(script,
                                  job_name="run_datruf",
                                  n_core=args.n_core,
                                  sync=False)
            elif args.job_scheduler == "slurm":
                script = slurm_nize(script,
                                    job_name="run_datruf",
                                    n_core=args.n_core,
                                    mem_per_cpu=40000,
                                    wait=False)
            f.write(script)

        run_command(f"{args.submit_job} {script_fname}")

    # Prepare a script for finalization of the task
    with open("finalize_datruf.sh", 'w') as f:
        f.write('\n'.join([f"cat {args.out_units_fname}.* > {args.out_units_fname}.cat",
                           f"cat {args.out_main_fname}.* > {args.out_main_fname}.cat",
                           f"awk -F'\\t' 'NR == 1 {{print $0}} $1 != \"\" {{print $0}}' {args.out_main_fname}.cat > {args.out_main_fname}",
                           f"rm {args.out_units_fname}.*; rm {args.out_main_fname}.*"]) + '\n')

    logger.info("Run `$ bash finalize_datruf.sh` after finishing all jobs")


def load_args():
    parser = argparse.ArgumentParser(
        description="Distribute datruf_run.py jobs using a job scheduler.")

    parser.add_argument(
        "db_file",
        help="DAZZ_DB file")

    parser.add_argument(
        "las_file",
        help="output of TANmask of the modified DAMASTER package")

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
        "-m",
        "--out_main_fname",
        type=str,
        default="datruf_result",
        help=("Write main results to this file. [datruf_result]"))

    parser.add_argument(
        "-u",
        "--out_units_fname",
        type=str,
        default="datruf_units",
        help=("Write unit sequences to this file. [datruf_units]"))

    parser.add_argument(
        "--only_interval",
        action="store_true",
        default=False,
        help=("Stop calculation just after obtaining TR intervals. Since "
              "filtering of TRs by CV of its unit lengths is not applied, "
              "(intervals of) TRs with short (<50 bp) units will be output, "
              "unlike the default mode. [False]"))

    parser.add_argument(
        "--on_the_fly",
        action="store_true",
        default=False,
        help=("Generate dump data for each read on the fly. This mode is very "
              "slow and used only when whole data are huge and you just want "
              "to look at results of only several reads. [False]"))

    parser.add_argument(
        "-D",
        "--debug_mode",
        action="store_true",
        default=False,
        help=("Run in debug mode. [False]"))

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

    return parser.parse_args()


if __name__ == "__main__":
    main()
