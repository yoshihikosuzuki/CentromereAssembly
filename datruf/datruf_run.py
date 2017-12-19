import argparse
import os
from collections import defaultdict
import numpy as np
import pandas as pd
from multiprocessing import Pool

from datruf_io import (load_dbdump,
                       load_ladump,
                       load_tr_intervals,
                       load_alignments,
                       load_paths)

from datruf_core import (calc_cover_set,
                         calc_min_cover_set)

from datruf_utils import (run_command,
                          add_element)


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

    columns = ("dbid", "start", "end", "unit length")
    columns_only_interval = ("dbid", "start", "end")

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

        result = defaultdict(dict)
        count = 0
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

            if self.only_interval:
                intervals = sorted([(x[2], x[1])
                                    for x in list(self.min_cover_set)])
                for start, end in intervals:
                    add_element(result,
                                Runner.columns_only_interval,
                                count,
                                (self.read_id,
                                 start,
                                 end))
                    count += 1
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

                # Add result into pd.DataFrame
                add_element(result, Runner.columns, count,
                            (self.read_id,
                             path.bb,
                             path.ae,
                             path.unit_len["mean"]))   # NOTE: border, mean, median
                count += 1

        result = pd.DataFrame.from_dict(result)
        return result


def run_runner(r):
    return r.run()


def main():
    args = load_args()
    n_core = args.n_core
    out_file = args.out_file
    del args.n_core
    del args.out_file

    if n_core == 1:
        r = Runner(args)
        results = r.run()
    else:
        r_tmp = Runner(args)   # only for adjusting dbids
        start_dbid = r_tmp.start_dbid
        end_dbid = r_tmp.end_dbid
        del r_tmp

        if n_core > end_dbid - start_dbid + 1:
            n_core = end_dbid - start_dbid + 1

        # Generate split Runners
        unit_nread = -(-(end_dbid - start_dbid + 1) // n_core)
        r = []
        for i in range(n_core):
            args.start_dbid = start_dbid + i * unit_nread
            args.end_dbid = min(end_dbid,
                                start_dbid + (i + 1) * unit_nread - 1)
            r.append(Runner(args))

        # Parallely execute and then merge results
        exe_pool = Pool(n_core)
        results = pd.DataFrame()
        for result in exe_pool.imap(run_runner, r):
            results = pd.concat([results, result])

        results = (results.sort_values(by="dbid", kind="mergesort")
                   .reset_index(drop=True))

    results.loc[:, Runner.columns].to_csv(out_file, sep="\t")


def load_args():
    parser = argparse.ArgumentParser(description="Run datruf with many reads")

    parser.add_argument("db_file",
                        help="DAZZ_DB file")

    parser.add_argument("las_file",
                        help="output of TANmask")

    parser.add_argument("--dbdump",
                        type=str,
                        default="datander_dbdump",
                        help=("Output of `DBdump -r -h -mtan <db_file>`. "
                              "This will be automatically generated if it "
                              "does not exist. In <on_the_fly> mode, this "
                              "is not used. [datander_dbdump]"))

    parser.add_argument("--ladump",
                        type=str,
                        default="datander_ladump",
                        help=("Output of `LAdump -c <db_file> <las_file>`. "
                              "This will be automatically generated if it "
                              "does not exist. In <on_the_fly> mode, this "
                              "is not used. [datander_ladump]"))

    parser.add_argument("--out_file",
                        type=str,
                        default="datruf_result",
                        help="A file for outputting results. [datruf_result]")

    parser.add_argument("--start_dbid",
                        type=int,
                        default=1,
                        help=("Start read ID, which is used in DAZZ_DB. "
                              "Set <= 1 for starting from the first read. "
                              "[1]"))

    parser.add_argument("--end_dbid",
                        type=int,
                        default=-1,
                        help=("End read ID. Set < 1 for ending at the last "
                              "read. [-1]"))

    parser.add_argument("--n_core",
                        type=int,
                        default=1,
                        help=("Degree of parallelization. [1]"))

    parser.add_argument("--only_interval",
                        action="store_true",
                        default=False,
                        help=("Stop calculation just after obtaining TR "
                              "intervals. [False]"))

    parser.add_argument("--on_the_fly",
                        action="store_true",
                        default=False,
                        help=("Generate dump data for each read on the fly. "
                              "This mode is very slow and used only when "
                              "whole data are huge and one just wants to "
                              "look at results of only several reads. "
                              "[False]"))

    return parser.parse_args()


if __name__ == "__main__":
    main()
