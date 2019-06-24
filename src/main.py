import argparse
from dataclasses import dataclass, field, InitVar
from os.path import isfile
import toml
from BITS.util.scheduler import Scheduler
from BITS.util.log import print_log
from .datander import run_datander
from .datruf import run_datruf

tasks = ['datander', 'datruf', 'dacmaster', 'dalayout']


@dataclass(repr=False, eq=False)
class VCA:
    '''Config file must be TOML-formatted.'''
    config_fname : InitVar[str]
    config       : dict         = field(init=False)
    scheduler    : Scheduler    = field(init=False, default=None)

    def __post_init__(self, config_fname):
        assert isfile(config_fname), f'No config file: {config_fname}'
        self.config = toml.load(config_fname)
        assert 'db_prefix' in self.config, 'Config must have a value "db_prefix"'

        if self.config.get('job_scheduler', {}).get('enabled', False):
            self.scheduler = Scheduler(**self.config['job_scheduler'].get('params', {}))

    def run(self):
        for task in tasks:
            self._run_datander()

    @print_log("datander")
    def _run_datander(self):
        # TODO: check skip condition
        run_datander(self.config['db_prefix'],
                     **self.config.get('datander', {}),
                     scheduler=self.scheduler)

    @print_log("datruf")
    def _run_datruf(self):
        run_datruf(f'{self.config["db_prefix"]}.db',
                   f'TAN.{self.config["db_prefix"]}.las',
                   **self.config.get('datruf', {}),
                   scheduler=self.scheduler)

    @print_log("dacmaster")
    def _run_dacmaster(self):
        pass

    @print_log("dalayout")
    def _run_dalayout(self):
        pass


def load_args():
    p = argparse.ArgumentParser(description="VCA: Vertebrate Centromere Assembler.")
    p.add_argument("config_fname", type=str, default="config",
                   help="TOML-formatted config file. [config]")
    return p.parse_args()


def main():
    args = load_args()
    vca = VCA(args.config_fname)
    vca.run()
