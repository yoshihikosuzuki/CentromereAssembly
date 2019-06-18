import argparse
from dataclasses import dataclass, field
from typing import Any
from os.path import isfile
from logzero import logger
from BITS.util.io import load_config
from BITS.util.proc import run_command
from BITS.util.scheduler import Scheduler
from BITS.util.log import print_log
from .datander import run_datander

tasks = ["all", "datander", "datruf", "dacmaster", "dalayout"]


def prepare_scheduler(config):
    if not config.getboolean('JOB_SCHEDULER', 'USE_SCHEDULER'):
        return None
    config_js = config['JOB_SCHEDULER']
    return Scheduler(config_js['SCHEDULER_NAME'],
                     config_js['SUBMIT_COMMAND'],
                     config_js['QUEUE_NAME'] if 'QUEUE_NAME' in config_js else None)


@dataclass(repr=False, eq=False)
class VCA:
    task_name    : str       = "all"
    config_fname : str       = "config"
    config       : Any       = field(init=False)
    scheduler    : Scheduler = field(init=False)

    def __post_init__(self):
        self.config = load_config(self.config_fname)
        self.scheduler = prepare_scheduler(self.config)

    def run(self):
        if self.task_name in ["all", "datander"]:
            self._run_datander()
        if self.task_name in ["all", "datruf"]:
            self._run_datruf()
        if self.task_name in ["all", "dacmaster"]:
            self._run_dacmaster()
        if self.task_name in ["all", "dalayout"]:
            self._run_dalayout()

    @print_log("datander")
    def _run_datander(self):
        run_datander(self.config['DAZZ_DB']['DB_PREFIX'],
                     int(self.config['DATANDER']['N_CORE']),
                     self.scheduler)

    @print_log("datruf")
    def _run_datruf(self):
        run_command(' '.join([f"datruf.py",
                              f"-n {config.get('DATRUF', 'N_CORE')}",
                              f"-p {config.get('DATRUF', 'N_DISTRIBUTE')}" if use_scheduler(config) else '',
                              f"{config.get('DAZZ_DB', 'DB_PREFIX')}.db",
                              f"TAN.{config.get('DAZZ_DB', 'DB_PREFIX')}.las"]))

    @print_log("dacmaster")
    def _run_dacmaster(self):
        run_command(' '.join([f"dacmaster.py",
                              f"-n {config.get('DATANDER', 'N_CORE')}",
                              scheduler_options,
                              f"{config.get('DAZZ_DB', 'DB_PREFIX')}"]))

    @print_log("dalayout")
    def _run_dalayout(self):
        run_command(' '.join([f"dalayout.py",
                              f"-n {config.get('DATANDER', 'N_CORE')}",
                              scheduler_options,
                              f"{config.get('DAZZ_DB', 'DB_PREFIX')}"]))


def load_args():
    p = argparse.ArgumentParser(description="VCA: Vertebrate Centromere Assembler.")

    p.add_argument("task_name", nargs='?', type=str, default="all",
                   help=f"Task name. {{{', '.join(tasks)}}} [all]")

    p.add_argument("-c", "--config_fname", type=str, default="config",
                   help="Config file name. [config]")

    args = p.parse_args()
    assert args.task_name in set(tasks), f"Invalid task name: {args.task_name}"
    assert isfile(args.config_fname), f"No config file: {args.config_fname}"
    logger.info(f"Starting task {args.task_name} with config file {args.config_fname}")
    return args


def main():
    args = load_args()
    vca = VCA(args.task_name, args.config_fname)
    vca.run()
