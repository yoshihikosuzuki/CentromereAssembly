import argparse
import os
import numpy as np
import pandas as pd

from datruf_io import (load_dbdump,
                       load_ladump,
                       load_tr_intervals,
                       load_alignments,
                       load_paths)

from datruf_core import (calc_cover_set,
                         calc_min_cover_set)

from datruf_utils import (run_command)


class Runner():
    """
    First, load all dump data of reads specified by <start_dbid> and
    <end_dbid>. If <on_the_fly> option is set, this loading is skipped
    and dump data will be loaded by directly running DBdump and LAdump
    at the time they are needed. <on_the_fly> mode is very slow due to
    its IO processing, so one should use this only when one just wants
    to look at results of only several reads without generating dbdump
    and ladump files.

    Then, for each read ID from <start_dbid> to <end_dbid>, apply datruf
    algorithm. Alignment paths are loaded on-the-fly in both modes due
    to its huge size.
    """

    def __init__(self, args):
        for attr in list(vars(args).keys()):
            setattr(self, attr, getattr(args, attr))

        # Adjust the start/end read IDs
        self._get_max_dbid()
        if self.start_dbid < 1:
            self.start_dbid = 1
        if self.end_dbid < 1 or self.end_dbid > self.max_dbid:
            self.end_dbid = self.max_dbid

    def _get_max_dbid(self):
        command = ("DBdump %s | awk 'NR == 1 {print $3}'") % self.db_file
        self.max_dbid = int(run_command(command).strip())

    def _check_dump(self):
        if not os.path.isfile(self.dbdump):
            command = ("DBdump -r -h -mtan %s > %s"
                       % (self.db_file, self.dbdump))
            run_command(command)

        if not os.path.isfile(self.ladump):
            command = ("LAdump -c %s %s > %s"
                       % (self.db_file, self.las_file, self.ladump))
            run_command(command)

    def run(self):
        if not self.on_the_fly:   # load dump data beforehand
            self._check_dump()
            self.tr_intervals_all = load_dbdump(self)
            self.alignments_all = load_ladump(self)

        for read_id in range(self.start_dbid, self.end_dbid + 1):
            self.read_id = read_id

            #print(self.read_id)

            # [(start, end), ...]
            self.tr_intervals = load_tr_intervals(self)
            if len(self.tr_intervals) == 0:
                continue

            # pd.DataFrame([abpos, aepos, bbpos, bepos, distance])
            # sorted by distance -> abpos
            self.alignments = load_alignments(self)

            # set((bb, ae, ab, be), ...)
            self.cover_set = calc_cover_set(self)

            # set((ab, ae, bb, be), ...)
            self.min_cover_set = calc_min_cover_set(self.cover_set)
            if len(self.min_cover_set) == 0:
                continue

            # [Path, ...]
            self.paths = load_paths(self)   # only alignments are loaded

            #print("---")
            #print(self.read_id)
            #print(self.tr_intervals)
            #print(self.alignments)
            #print(self.cover_set)
            #print(self.min_cover_set)

            for path in self.paths:
                path.split_alignment()
                print(self.read_id, path.bb, path.ae, path.unit_len["border"], path.unit_len["mean"], path.unit_len["median"])   # TODO: use dataframe


def main():
    r = Runner(load_args())
    r.run()
    # TODO: parallelization (both inside this (multiprocessing) and outside this (by shell script))


def load_args():
    parser = argparse.ArgumentParser(description="Run datruf with many reads")

    parser.add_argument("db_file",
                        help="DAZZ_DB file")

    parser.add_argument("las_file",
                        help="output of TANmask")

    parser.add_argument("--dbdump",
                        type=str,
                        default="datander_dbdump",
                        help="output of `DBdump -r -h -mtan <db_file>`")

    parser.add_argument("--ladump",
                        type=str,
                        default="datander_ladump",
                        help="output of `LAdump -c <db_file> <las_file>`")

    parser.add_argument("--out_file",
                        type=str,
                        default="datander_result",
                        help="")

    parser.add_argument("--start_dbid",
                        type=int,
                        default=18,
                        help=("Start read ID, which is used in DAZZ_DB."
                              "<= 1 for starting from the first read. [1]"))

    parser.add_argument("--end_dbid",
                        type=int,
                        default=30,
                        help=("End read ID. < 1 for ending at the last read. "
                              "[-1]"))

    parser.add_argument("--on_the_fly",
                        action="store_true",
                        default=False,
                        help=("Generate dump data for each read on the fly."
                              "This mode is very slow and used only when "
                              "whole data are huge and one just wants to "
                              "look at results of only several reads. "
                              "[False]"))

    return parser.parse_args()


if __name__ == "__main__":
    main()
