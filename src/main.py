import argparse
import configparser
from logzero import logger
from BITS.utils import run_command, print_log


@print_log("datander")
def run_datander(config, job_scheduler_options):
    run_command(' '.join([f"datander.py",
                          f"-n {config.get('DATANDER', 'N_CORE')}",
                          job_scheduler_options,
                          f"{config.get('DAZZ_DB', 'DB_PREFIX')}"]))


@print_log("datruf")
def run_datruf(config, job_scheduler_options):
    run_command(' '.join([f"datruf.py",
                          f"-n {config.get('DATANDER', 'N_CORE')}",
                          job_scheduler_options,
                          f"{config.get('DAZZ_DB', 'DB_PREFIX')}"]))


@print_log("dacmaster")
def run_dacmaster(config, job_scheduler_options):
    run_command(' '.join([f"dacmaster.py",
                          f"-n {config.get('DATANDER', 'N_CORE')}",
                          job_scheduler_options,
                          f"{config.get('DAZZ_DB', 'DB_PREFIX')}"]))


@print_log("dalayout")
def run_dalayout(config, job_scheduler_options):
    run_command(' '.join([f"dalayout.py",
                          f"-n {config.get('DATANDER', 'N_CORE')}",
                          job_scheduler_options,
                          f"{config.get('DAZZ_DB', 'DB_PREFIX')}"]))


def main():
    args = load_args()
    config = configparser.ConfigParser()
    config.read(args.config_fname)

    # Prepare options for job scheduler
    job_scheduler_options = ""
    if 'USE_SCHEDULER' in config['JOB_SCHEDULER'] and config.getboolean('JOB_SCHEDULER', 'USE_SCHEDULER'):
        job_scheduler_options = [f"--job_scheduler {config.get('JOB_SCHEDULER', 'SCHEDULER_NAME')}",
                                 f"--submit_command {config.get('JOB_SCHEDULER', 'SUBMIT_COMMAND')}"]
        if 'QUEUE_NAME' in config['JOB_SCHEDULER']:
            job_scheduler_options.append(f"--queue_name {config.get('JOB_SCHEDULER', 'QUEUE_NAME')}")
        if 'TIME_LIMIT' in config['JOB_SCHEDULER']:
            job_scheduler_options.append(f"--time_limit {config.get('JOB_SCHEDULER', 'TIME_LIMIT')}")
        if 'MEM_LIMIT' in config['JOB_SCHEDULER']:
            job_scheduler_options.append(f"--mem_limit {config.get('JOB_SCHEDULER', 'MEM_LIMIT')}")
        job_scheduler_options = ' '.join(job_scheduler_options)

    # Run submodule(s)
    if args.task_name in ("datander", "all"):
        run_datander(config, job_scheduler_options)
    if args.task_name in ("datruf", "all"):
        run_datruf(config, job_scheduler_options)
    if args.task_name in ("dacmaster", "all"):
        run_dacmaster(config, job_scheduler_options)
    if args.task_name in ("dalayout", "all"):
        run_dalayout(config, job_scheduler_options)


def load_args():
    p = argparse.ArgumentParser(description="Vertebrate Centromere Assembler.")

    p.add_argument("task_name",
                   nargs='?',
                   type=str,
                   default="all",
                   help="Task name. This must be one of {'all', 'datruf', 'dacmaster', 'dalayout'}. [all]")

    p.add_argument("-c",
                   "--config_fname",
                   type=str,
                   default="config",
                   help="Config file name. [config]")

    args = p.parse_args()
    assert args.task_name in set(["all", "datander", "datruf", "dacmaster", "dalayout"]), "Invalid task name"
    return args
