import numpy as np
from logzero import logger
from BITS.utils import run_command
from BITS.scheduler import Scheduler
from .run import Runner


def main():
    args = load_args()

    # Use Runner instance just for some preparation of data
    r = Runner(args)   # args are now variables of <r>
    r._check_dump()

    n_dbid_part = -(-r.end_dbid // args.n_distribute)
    n_digit = int(np.log10(args.n_distribute) + 1)

    # Submit scripts each of which runs datruf_run.py
    logger.info("Scattering jobs. Intermediate files are stored in datruf/")
    run_command("mkdir -p datruf; rm -f datruf/*")
    s = Scheduler(args.job_scheduler,
                  args.submit_command,
                  args.queue_name,
                  "datruf/log")
    jids = []
    for i in range(args.n_distribute):
        index = str(i + 1).zfill(n_digit)
        jids.append(s.submit(' '.join([f"datruf_run.py",
                                       f"-s {r.start_dbid + i * n_dbid_part}",
                                       f"-e {r.start_dbid + (i + 1) * n_dbid_part - 1}",
                                       f"-m datruf/{args.out_main_fname}.{index}",
                                       f"-u datruf/{args.out_units_fname}.{index}",
                                       f"{'--only_interval' if args.only_interval else ''}",
                                       f"{'-D' if args.debug_mode else ''}",
                                       f"-n {args.n_core}",
                                       f"{args.db_file}",
                                       f"{args.las_file}"]),
                             f"datruf/run_datruf.{args.job_scheduler}.{index}",
                             job_name="datruf_dist",
                             n_core=args.n_core))

    # Merge the results
    logger.info("Waiting for all jobs to be finished...")
    s.submit('\n'.join([f"cat datruf/{args.out_main_fname}.* > datruf/{args.out_main_fname}.cat",
                        f"cat datruf/{args.out_units_fname}.* > datruf/{args.out_units_fname}.cat",
                        f"awk -F'\\t' 'BEGIN {{count = 0}} NR == 1 {{print $0}} $1 != \"\" {{printf count; for (i = 2; i <= NF; i++) {{printf \"\\t\" $i}}; print \"\"; count++}}' datruf/{args.out_main_fname}.cat > {args.out_main_fname}",
                        f"awk -F'\\t' 'BEGIN {{count = 0}} NR == 1 {{print $0}} $1 != \"\" {{printf count; for (i = 2; i <= NF; i++) {{printf \"\\t\" $i}}; print \"\"; count++}}' datruf/{args.out_units_fname}.cat > {args.out_units_fname}"]),
             f"datruf/finalize_datruf.{args.job_scheduler}",
             job_name="datruf_finalize",
             depend=jids,
             wait=True)


def load_args():
    import argparse
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
        "-n",
        "--n_core",
        type=int,
        default=1,
        help=("Degree of parallelization in each distributed job. [1]"))

    parser.add_argument(
        "-D",
        "--debug_mode",
        action="store_true",
        default=False,
        help=("Run datruf in debug mode. [False]"))

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
        "--queue_name",
        type=str,
        default=None,
        help="Name of queue (SGE) or partition (SLURM) to which jobs are submitted. [None]")

    return parser.parse_args()


if __name__ == "__main__":
    main()
