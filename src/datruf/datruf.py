from dataclasses import dataclass
import argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool
from interval import interval
from logzero import logger
from BITS.util.interval import intvl_len, subtract_intvl
from BITS.util.proc import run_command
from BITS.util.scheduler import Scheduler
from ..types import TRRead

dir_name     = "datruf"
script_scatter_fname = f"{dir_name}/scatter.sh"
script_gather_fname = f"{dir_name}/gather.sh"
out_fname    = f"{dir_name}/datruf_result"
log_fname    = f"{dir_name}/log"


@dataclass(eq=False)
class DatrufRunner:
    """Entry point of datruf, which detects units of TRs using the result of datander.

    Positional arguments:
      - db_fname  <str> : DAZZ_DB file
      - las_fname <str> : Output of datander. These files must be in CWD
    
    Optional arguments:
      - n_core       <int>       [1]               : Number of cores used in a single job of datrud
      - n_distribute <int>       [1]               : Number of jobs distributed in datruf
      - scheduler    <Scheduler> [None]            : Scheduler object
    """
    db_fname     : str
    las_fname    : str
    n_core       : int       = 1
    n_distribute : int       = 1
    scheduler    : Scheduler = None

    def __post_init__(self):
        run_command(f"mkdir -p {dir_name}; rm -f {dir_name}/*")

    def run(self):
        pass
