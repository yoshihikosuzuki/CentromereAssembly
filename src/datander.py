from os.path import join
from dataclasses import dataclass
from logzero import logger
from BITS.util.proc import run_command
from BITS.util.scheduler import Scheduler


@dataclass(eq=False)
class DatanderRunner:
    """Entry point of datander, a commandline tool for detecting tandem repeat regions from (noisy) reads.
    In VCA, slightly customized datander is used (alignment will NOT be extended to the ends of a read).

    Positional arguments:
      - db_prefix <str> : Prefix of the DB file created with DAZZ_DB. DB file must be in CWD

    Optional arguments:
      - n_core       <int>       [1]                 : Number of cores used in datader
      - scheduler    <Scheduler> [None]              : Scheduler object
      - dir_name     <str>       ["datander"]        : Directory name to which results will be output
      - script_fname <str>       ["run_datander.sh"] : File name of the script used for executing HPC.TANmask
      - log_fname    <str>       ["log"]             : Log file name
    """
    db_prefix    : str
    n_core       : int       = 1
    scheduler    : Scheduler = None
    dir_name     : str       = "datander"
    script_fname : str       = "run_datander.sh"
    log_fname    : str       = "log"

    def __post_init__(self):
        self.script_fname = join(self.dir_name, self.script_fname)
        self.log_fname = join(self.dir_name, self.log_fname)

        run_command(f"rm -f .{self.db_prefix}.*.tan.* .{self.db_prefix}.tan.* TAN.*")
        run_command(f"mkdir -p {self.dir_name}; rm -f {self.dir_name}/*")

    def run(self):
        # Prepare a script to run datander
        script = run_command(f"HPC.TANmask -T{self.n_core} {self.db_prefix}.db")
        if calc_n_blocks(f"{self.db_prefix}.db") > 1:
            script += '\n'.join([f"Catrack -v {self.db_prefix} tan",
                                 f"rm .{self.db_prefix}.*.tan.*"])

        # Run the script
        if self.scheduler is None:
            with open(self.script_fname, "w") as f:
                f.write(f"{script}\n")
                run_command(f"bash {self.script_fname} > {self.log_fname} 2>&1")
        else:
            self.scheduler.submit(script,
                                  self.script_fname,
                                  job_name=self.dir_name,
                                  log_fname=self.log_fname,
                                  n_core=self.n_core,
                                  wait=True)


def calc_n_blocks(db_fname):
    """Extract the number of blocks from the db file."""
    with open(db_fname, 'r') as f:
        for line in f:
            if line.startswith("blocks"):
                return int(line.split('=')[1].strip())
    logger.error(f"No information on the number of blocks in {db_fname}")
