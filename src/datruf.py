import argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool
from interval import interval
from logzero import logger
from BITS.scheduler import Scheduler
from BITS.interval import interval_len, subtract_interval
from BITS.utils import run_command, debug_mode
from .scheduler_args import add_scheduler_args

dir_name = "datruf"
out_fname = "datruf_result"


def filter_alignments(tr_intervals, alignments, min_len=1000):
    # NOTE: <min_len> defines the required overlap length with yet uncovered TR region
    uncovered = interval(*tr_intervals)
    inners = set()
    for ab, ae, bb, be, distance, slope in alignments:   # in order of distance
        if interval_len(uncovered) < min_len:
            break
        intersect = uncovered & interval[bb, ae]
        uncovered = subtract_interval(uncovered, interval[bb, ae])
        if (interval_len(intersect) >= min_len
            and 0.95 <= slope <= 1.05   # eliminate abnornal slope
            and ab <= be):   # at least duplication
            inners.add((ab, ae, bb, be))   # TODO: add only intersection is better?
    logger.debug(f"inners: {inners}")
    return inners


def load_paths(args, read_id, alignment_set):
    # NOTE: Since alignment path information is very large, load for each read on demand

    def find_boundary(aseq, bseq):
        # NOTE: '[' and ']' are alignment boundary, '...' is read boundary
        start = aseq.find('[') + 1
        if start == 0:
            while aseq[start] != '.' and bseq[start] != '.':
                start += 1
            while aseq[start] == '.' or bseq[start] == '.':
                start += 1
        end = aseq.rfind(']')
        if end == -1:
            end = len(aseq)
            while aseq[end - 1] != '.' and bseq[end - 1] != '.':
                end -= 1
            while aseq[end - 1] == '.' or bseq[end - 1] == '.':
                end -= 1
        return start, end

    def convert_symbol(aseq, bseq, symbol):
        return ''.join(['M' if c == '|'
                        else 'I' if aseq[i] == '-'
                        else 'D' if bseq[i] == '-'
                        else 'X'
                        for i, c in enumerate(symbol)])

    if len(alignment_set) == 0:
        return []

    out = run_command(f"LAshow4pathplot -a {args.db_file} {args.las_file} {read_id} | sed 's/,//g' | awk 'BEGIN {{first = 1}} NF == 7 {{if (first == 1) {{first = 0}} else {{printf(\"%s\\n%s\\n%s\\n%s\\n\", header, aseq, bseq, symbol)}}; header = $0; aseq = \"\"; bseq = \"\"; symbol = \"\"; count = 0;}} NF < 7 {{if (count == 0) {{aseq = aseq $0}} else if (count == 1) {{symbol = symbol $0}} else {{bseq = bseq $0}}; count++; count %= 3;}} END {{printf(\"%s\\n%s\\n%s\\n%s\\n\", header, aseq, bseq, symbol)}}'").strip().split('\n')

    ret = []   # [(ab, ae, bb, be, FlattenCigar), ...]
    for header, aseq, bseq, symbol in zip(*([iter(out)] * 4)):    # split every 4 lines (= single entry)
        alignment = tuple(map(int, header.replace(' ', '').split('\t')[2:6]))   # ab, ae, bb, be
        if alignment in alignment_set:
            ret.append((*alignment, convert_symbol(*map(lambda x: x[slice(*find_boundary(aseq, bseq))],
                                                        [aseq, bseq, symbol]))))
    logger.debug(f"read {read_id}: {len(ret)} paths loaded")
    return sorted(ret, key=lambda x: x[2])   # sort by bb


def split_tr(ab, bb, fcigar):
    """
    Split TR interval [bb, ae] into unit intervals.
    """

    apos, bpos = ab, bb
    unit_intvls = [(bpos, apos)]
    # Iteratively find max{ax} such that bx == (last unit end)
    for i, c in enumerate(fcigar):
        if c != 'I':
            apos += 1
        if c != 'D':
            bpos += 1
        if bpos == unit_intvls[-1][1] and (i == len(fcigar) - 1 or fcigar[i + 1] != 'D'):
            unit_intvls.append((unit_intvls[-1][1], apos))
    return unit_intvls


def find_units_single(args, read_id, tr_intervals, alignments, max_cv=0.1):
    ret = []
    for ab, ae, bb, be, fcigar in load_paths(args, read_id, filter_alignments(tr_intervals, alignments)):
        unit_intvls = split_tr(ab, bb, fcigar)
        if len(unit_intvls) == 1:   # at least duplication
            continue
        unit_lens = [t - s for s, t in unit_intvls]
        cv_len = round(np.std(unit_lens, ddof=1) / np.mean(unit_lens), 3)
        if cv_len >= max_cv:   # probably due to STR
            continue
        mean_len = int(np.mean(unit_lens))
        ret.append((read_id, bb, ae, mean_len, cv_len, unit_intvls))
    return ret


def find_units_mult(args, read_ids, tr_intervals_all, alignments_all):
    return [find_units_single(args, read_id, tr_intervals_all[read_id], alignments_all[read_id])
            for read_id in read_ids]


def load_dumps(args):
    # Extract data from DBdump's and LAdump's output
    return ({read_id: list(df.apply(lambda d: tuple(d[["start", "end"]]), axis=1))
             for read_id, df
             in (pd.DataFrame([list(map(int, x.split('\t')))
                               for x in run_command(f"DBdump -rh -mtan {args.db_file} {args.start_dbid}-{args.end_dbid} | awk '$1 == \"R\" {{dbid = $2}} $1 == \"T0\" && $2 > 0 {{for (i = 1; i <= $2; i++) printf(\"%s\\t%s\\t%s\\n\", dbid, $(2 * i + 1), $(2 * i + 2))}}'").strip().split('\n')],
                              columns=("dbid", "start", "end")) \
                 .groupby("dbid"))},
            {read_id: list(df.sort_values(by="abpos", kind="mergesort") \
                           .sort_values(by="distance", kind="mergesort") \
                           .apply(lambda d: tuple((*map(int, d[["abpos", "aepos", "bbpos", "bepos", "distance"]]), d["slope"])), axis=1))
             for read_id, df
             in (pd.DataFrame([list(map(int, x.split('\t')))
                               for x in run_command(f"LAdump -c {args.db_file} {args.las_file} {args.start_dbid}-{args.end_dbid} | awk '$1 == \"P\" {{dbid = $2}} $1 == \"C\" {{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\", dbid, $2, $3, $4, $5)}}'").strip().split('\n')],
                              columns=("dbid", "abpos", "aepos", "bbpos", "bepos")) \
                 .assign(distance=lambda x: x["abpos"] - x["bbpos"]) \
                 .assign(slope=lambda x: ((x["aepos"] - x["abpos"]) / (x["bepos"] - x["bbpos"])).round(3)) \
                 .groupby("dbid"))})


def find_units(args):
    tr_intervals_all, alignments_all = load_dumps(args)
    read_ids = sorted(tr_intervals_all.keys())

    results = {}
    index = 0
    if args.n_core == 1:
        logger.debug("single core")
        ret = find_units_mult(args, read_ids, tr_intervals_all, alignments_all)
        for r in ret:
            if len(r) != 0:
                results[index] = r
                index += 1
    else:
        logger.debug("multi core")
        # split list
        n_per_core = -(-len(read_ids) // args.n_core)
        list_args = [(args,
                      read_ids[i * n_per_core : (i + 1) * n_per_core],
                      tr_intervals_all,
                      alignments_all)
                     for i in range(args.n_core)]
        with Pool(args.n_core) as pool:
            for rets in pool.starmap(find_units_mult, list_args):
                for ret in rets:
                    for r in ret:
                        if len(r) != 0:
                            results[index] = r
                            index += 1

    pd.DataFrame.from_dict(results,
                           orient="index",
                           columns=("read_id",
                                    "start",
                                    "end",
                                    "mean_ulen",
                                    "cv_ulen",
                                    "units")) \
                .to_csv(out_fname if args.index is None else f"{dir_name}/{out_fname}.{args.index}", sep='\t')


def concat_df(dir_name, prefix, sep='\t', index_col=0):
    return pd.concat([pd.read_csv(fname, sep='\t', index_col=index_col)
                      for fname in run_command(f"find {dir_name} -name '{prefix}.*' | sort").strip().split('\n')]) \
             .reset_index(drop=True)


def main():
    args = load_args()
    if args.end_dbid <= 0:   # Set <end_dbid> as the last read if not specified
        args.end_dbid = int(run_command(f"DBdump {args.db_file} | awk 'NR == 1 {{print $3}}'").strip())
    run_command(f"mkdir -p {dir_name}; rm -f {dir_name}/*")

    if args.job_scheduler is None:
        logger.debug("Without job scheduler")
        find_units(args)
    else:
        logger.debug("With job scheduler")
        # Split tasks and submit jobs
        n_part = -(-(args.end_dbid - args.start_dbid + 1) // args.n_distribute)
        indices = [str(i + 1).zfill(int(np.log10(args.n_distribute) + 1))
                   for i in range(args.n_distribute)]
        s = Scheduler(args.job_scheduler,
                      args.submit_command,
                      args.queue_name,
                      f"{dir_name}/log")
        jids = [s.submit(' '.join([f"python -m dacembler.datruf",
                                   f"-s {args.start_dbid + i * n_part}",
                                   f"-e {min([args.start_dbid + (i + 1) * n_part - 1, args.end_dbid])}",
                                   f"-n {args.n_core}",
                                   f"{'-D' if args.debug_mode else ''}",
                                   f"--index {index}",
                                   f"{args.db_file}",
                                   f"{args.las_file}"]),
                         f"{dir_name}/run.sh.{index}",
                         job_name="datruf_distribute",
                         n_core=args.n_core)
                for i, index in enumerate(indices)]

        # Merge the results
        logger.info("Waiting for all distributed jobs to be finished...")
        s.submit("sleep 1s",
                 f"{dir_name}/gather.sh",
                 job_name="datruf_gather",
                 depend=jids,
                 wait=True)
        concat_df(dir_name, out_fname).round(3).to_csv(out_fname, sep='\t')


def main_distribute():
    find_units(load_args())


def load_args():
    p = argparse.ArgumentParser(description="Find tandem repeat units from PacBio reads.")

    p.add_argument("db_file",
                   help="DAZZ_DB file.")

    p.add_argument("las_file",
                   help="Output of TANmask of the modified DAMASTER package.")

    p.add_argument("-s",
                   "--start_dbid",
                   type=int,
                   default=1,
                   help="Read ID of DAZZ_DB from which datruf's computation starts. [1]")

    p.add_argument("-e",
                   "--end_dbid",
                   type=int,
                   default=-1,
                   help="Read ID of DAZZ_DB at which datruf's computation ends. -1 means the last read. [-1]")

    p.add_argument("-n",
                   "--n_core",
                   type=int,
                   default=1,
                   help="Degree of parallelization. [1]")

    p.add_argument("-p",
                   "--n_distribute",
                   type=int,
                   default=1,
                   help="Degree of parallelization. [1]")

    p.add_argument("-D",
                   "--debug_mode",
                   action="store_true",
                   default=False,
                   help="Show debug messages. [False]")

    p.add_argument("--index",
                   type=str,
                   default=None,
                   help="Internally used for distributed computation.")

    add_scheduler_args(p)
    args = p.parse_args()
    debug_mode(args.debug_mode)
    if args.n_distribute > 1:
        assert args.job_scheduler is not None, "--job_scheduler must be specified"
        assert args.submit_command is not None, "--submit_command must be specified"
    return args


if __name__ == "__main__":
    main_distribute()
