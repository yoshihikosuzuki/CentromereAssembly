from typing import List
from dataclasses import dataclass, field, InitVar
import numpy as np
import pandas as pd
from interval import interval
from logzero import logger
from BITS.interval import interval_len, subtract_interval
from BITS.seq import load_fasta, revcomp, homopolymer_compression
from BITS.utils import run_command


@dataclass
class Read:
    dbid: int
    header: str
    sequence: str

    @property
    def length(self):
        return len(self.sequence)


def load_db(db_file):
    # NOTE: header must not contain tab
    return {dbid: Read(dbid, *line.split('\t'))
            for dbid, line
            in enumerate(run_command(f"DBshow {db_file} | awk 'BEGIN {{first = 1}} {{if (substr($0, 1, 1) == \">\") {{if (first == 1) {{first = 0}} else {{printf(\"%s\\t%s\\n\", header, seq)}}; header = substr($0, 2); seq = "";}} else {{seq = seq $0}}}} END {{printf(\"%s\\t%s\\n\", header, seq)}}'").strip().split('\n'), start=1)}


@dataclass
class ReadInterval:
    read: Read
    start: int
    end: int

    @property
    def interval(self):
        return interval[self.start, self.end]

    @property
    def length(self):
        return self.end - self.start

    @property
    def sequence(self):
        return self.read.sequence[self.start:self.end]


@dataclass
class Unit(ReadInterval):
    repr_id: int = field(init=False)


@dataclass
class Alignment:
    ab: int
    ae: int
    bb: int
    be: int
    fcigar: str


@dataclass
class TR(ReadInterval):
    alignment: Alignment


def main():
    args = load_args()
    reads = load_db(args.db_file)
    for read_id in read_ids:
        inners.add(TR(reads[read_id], ))
