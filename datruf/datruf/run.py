import os.path
import numpy as np
import pandas as pd
import logging
import logzero
from logzero import logger
from BITS.utils import run_command
from .io import load_dbdump, load_ladump, load_tr_intervals, load_alignments, load_paths
from .core import calc_cover_set, calc_min_cover_set


class Runner():
    """
    First, load all dump data of reads specified by <start_dbid> and
    <end_dbid>. If <on_the_fly> option is set (not recommended in most cases),
    this loading is skipped and dump data will be loaded by executing DBdump
    and LAdump on demand.

    Then, for each read ID from <start_dbid> to <end_dbid>, apply the datruf
    algorithm in parallel (if <n_core> > 1). If <only_interval> option is set,
    all tandem repeat INTERVALs including those with short (< 50 bp) units will
    be output. Otherwise (i.e. by defalut), intervals and UNITs of tandem
    repeats after filtering by coefficient of variation of unit lengths will be
    output. This filtering usually removes tandem repeats with short units.
    """

    def __init__(self, args):
        for attr in list(vars(args).keys()):
            setattr(self, attr, getattr(args, attr))

        # handle with improper start/end read IDs specified
        self._get_max_dbid()
        if self.start_dbid < 1:
            self.start_dbid = 1
        if self.end_dbid < 1 or self.end_dbid > self.max_dbid:
            self.end_dbid = self.max_dbid

    def _get_max_dbid(self):
        command = f"DBdump {self.db_file} | awk 'NR == 1 {{print $3}}'"
        self.max_dbid = int(run_command(command).strip())

    def _check_dump(self):
        if not os.path.isfile(self.dbdump):
            logger.info("DBdump does not exist. Generating...")
            run_command(f"DBdump -r -h -mtan {self.db_file} > {self.dbdump}")
        if not os.path.isfile(self.ladump):
            logger.info("LAdump does not exist. Generating...")
            run_command(f"LAdump -c {self.db_file} {self.las_file} > {self.ladump}")

    def run(self):
        if not self.on_the_fly:
            # Generate dump data if not exist
            self._check_dump()

            # Load dump data
            self.tr_intervals_all = load_dbdump(self)
            if len(self.tr_intervals_all) == 0:
                return None
            self.alignments_all = load_ladump(self)

        trs, units = {}, {}   # to be output as <args.out_main_fname> and <args.out_units_fname>
        tr_index = unit_index = 0
        for read_id in range(self.start_dbid, self.end_dbid + 1):
            self.read_id = read_id

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
                intervals = sorted([(x[2], x[1])   # (bb, ae)
                                    for x in list(self.min_cover_set)])
                for start, end in intervals:
                    trs[tr_index] = [self.read_id,
                                     start,
                                     end,
                                     np.nan,
                                     np.nan,
                                     np.nan]
                    tr_index += 1
                continue

            # else: Default mode below

            # [Path, ...]
            self.paths = load_paths(self)
            
            for path_id, path in enumerate(self.paths):
                path.split_alignment()
                mean_len, cv = path.calc_unit_stat()
                trs[tr_index] = [self.read_id,
                                 path.bb,
                                 path.ae,
                                 len(path.unit_seqs),
                                 mean_len,
                                 cv]
                tr_index += 1

                # Also output unit sequences if the alignment is stable enough
                if np.isnan(cv) or cv >= 0.1:
                    continue
                for unit_id, unit_seq in enumerate(path.unit_seqs):
                    start, end = path.reflection_points[unit_id:unit_id + 2]
                    units[unit_index] = [self.read_id,
                                         path_id,
                                         unit_id,
                                         start,
                                         end,
                                         end - start,
                                         unit_seq]   # TODO: maybe no sequence is better (then load on demand)
                    unit_index += 1

        trs = pd.DataFrame.from_dict(trs,
                                     orient="index",
                                     columns=("read_id",
                                              "start",
                                              "end",
                                              "unit_count",
                                              "mean_unit_length",
                                              "unit_cv")).round(3)
        units = pd.DataFrame.from_dict(units,
                                       orient="index",
                                       columns=("read_id",
                                                "path_id",
                                                "unit_id",
                                                "start",
                                                "end",
                                                "length",
                                                "sequence"))
        return (trs, units)


def run_runner(r):
    return r.run()


def main():
    from multiprocessing import Pool

    args = load_args()
    n_core, out_main_fname, out_units_fname = args.n_core, args.out_main_fname, args.out_units_fname
    del args.n_core
    del args.out_main_fname
    del args.out_units_fname

    if n_core == 1:
        r = Runner(args)
        ret = r.run()
        if ret is None:
            return
        trs_all, units_all = ret
    else:
        # First of all, generate *_dump files if not exist
        r_tmp = Runner(args)
        r_tmp._check_dump()
        # Also obtain start/end read IDs in the whole data
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
        trs_all, units_all = pd.DataFrame(), pd.DataFrame()
        for ret in exe_pool.imap(run_runner, r):
            if ret is not None:
                trs, units = ret
                trs_all = pd.concat([trs_all, trs])
                units_all = pd.concat([units_all, units])
        exe_pool.close()

        if trs_all.shape[0] == 0:
            logger.info("No TRs. Exit without any output")
            return

        # Sort the results in the order of read_id
        trs_all = (trs_all.sort_values(by="read_id", kind="mergesort")
                   .reset_index(drop=True))

        if units_all.shape[0] != 0:
            units_all = (units_all.sort_values(by="read_id", kind="mergesort")
                         .reset_index(drop=True))

    trs_all.to_csv(out_main_fname, sep="\t")
    units_all.to_csv(out_units_fname, sep="\t")


def load_args():
    import argparse
    parser = argparse.ArgumentParser(
        description="Detect tandem repeat intervals and their unit sequences in PacBio reads.")

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
              "automatically generated if it does not exist. In "
              "<on_the_fly> mode, this is not used. [datander_dbdump]"))

    parser.add_argument(
        "-l",
        "--ladump",
        type=str,
        default="datander_ladump",
        help=("Output of `LAdump -c <db_file> <las_file>`. This will be "
              "automatically generated if it does not exist. In "
              "<on_the_fly> mode, this is not used. [datander_ladump]"))

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
        help=("Degree of parallelization. [1]"))

    parser.add_argument(
        "-D",
        "--debug_mode",
        action="store_true",
        default=False,
        help=("Run in debug mode. [False]"))

    args = parser.parse_args()
    if args.debug_mode:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    del args.debug_mode

    return args


if __name__ == "__main__":
    main()
