import argparse
from logzero import logger
from BITS.scheduler import Scheduler
from BITS.utils import run_command, debug_mode
from .scheduler_args import add_scheduler_args


def main():
    args = load_args()
    dir_name = "datander"

    # Prepare the directory
    run_command(f"rm -f .{args.db_prefix}.*.tan.* .{args.db_prefix}.tan.* TAN.*")
    run_command(f"mkdir -p {dir_name}; rm -f {dir_name}/*")

    # Run the script
    script = '\n'.join([run_command(f"HPC.TANmask -T{args.n_core} {args.db_prefix}.db"),
                        f"Catrack -v {args.db_prefix} tan",
                        f"rm .{args.db_prefix}.*.tan.*"])
    script_fname = f"{dir_name}/run_datander.sh"
    if args.job_scheduler is None:
        with open(script_fname, 'w') as f:
            f.write(f"{script}\n")
        run_command(f"bash {script_fname} > {dir_name}/log 2>&1")
    else:
        Scheduler(args.job_scheduler,
                  args.submit_command,
                  args.queue_name,
                  f"{dir_name}/log").submit(script,
                                            script_fname,
                                            job_name="run_datander",
                                            n_core=args.n_core,
                                            wait=True)


def load_args():
    p = argparse.ArgumentParser(description="Run datander.")

    p.add_argument("db_prefix",
                   help="Prefix name of a DAZZ_DB file.")

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

    add_scheduler_args(p)
    args = p.parse_args()
    debug_mode(args.debug_mode)
    return args
