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
        n_reads = int(run_command(f"DBdump {self.db_fname} | awk 'NR == 1 {{print $3}}'").strip())

        if self.scheduler is None:
            find_units(1, n_reads, self.db_fname, self.las_fname, self.n_core, out_fname)
        else:
            jids = []
            unit_n = -(-n_reads // self.n_distribute)
            for i in range(self.n_distribute):
                index = str(i + 1).zfill(int(np.log10(self.n_distribute) + 1))
                start = 1 + i * unit_n
                end = min([1 + (i + 1) * unit_n - 1, n_reads])
                script = (f"python -m vca.datruf {self.db_fname} {self.las_fname} {out_fname}.{index} "
                          f"{start} {end} {self.n_core}")
                
                jids.append(self.scheduler.submit(script,
                                                  f"{script_scatter_fname}.{index}",
                                                  job_name="datruf_scatter",
                                                  log_fname=log_fname,
                                                  n_core=self.n_core))
                
            # Merge the results
            logger.info("Waiting for all distributed jobs to be finished...")
            self.scheduler.submit("sleep 1s",
                                  script_gather_fname,
                                  job_name="datruf_gather",
                                  log_fname=log_fname,
                                  depend=jids,
                                  wait=True)
            concat_df(dir_name, out_fname).round(3).to_csv(out_fname, sep="\t")
            


def filter_alignments(tr_intervals, alignments, min_len=1000):
    # NOTE: <min_len> defines the required overlap length with yet uncovered TR region
    uncovered = interval(*tr_intervals)
    inners = set()
    for ab, ae, bb, be, distance, slope in alignments:   # in order of distance
        if intvl_len(uncovered) < min_len:
            break
        intersect = uncovered & interval[bb, ae]
        uncovered = subtract_intvl(uncovered, interval[bb, ae])
        if (intvl_len(intersect) >= min_len
            and 0.95 <= slope <= 1.05   # eliminate abnornal slope
            and ab <= be):   # at least duplication
            inners.add((ab, ae, bb, be))   # TODO: add only intersection is better?
    logger.debug(f"inners: {inners}")
    return inners


def load_paths(read_id, alignment_set, db_file, las_file):
    # NOTE: Since alignment path information is very large, load for each read on demand

    def find_boundary(aseq, bseq):
        # NOTE: "[" and "]" are alignment boundary, "..." is read boundary
        start = aseq.find("[") + 1
        if start == 0:
            while aseq[start] != "." and bseq[start] != ".":
                start += 1
            while aseq[start] == "." or bseq[start] == ".":
                start += 1
        end = aseq.rfind("]")
        if end == -1:
            end = len(aseq)
            while aseq[end - 1] != "." and bseq[end - 1] != ".":
                end -= 1
            while aseq[end - 1] == "." or bseq[end - 1] == ".":
                end -= 1
        return start, end

    def convert_symbol(aseq, bseq, symbol):
        return "".join(["M" if c == "|"
                        else "I" if aseq[i] == "-"
                        else "D" if bseq[i] == "-"
                        else "X"
                        for i, c in enumerate(symbol)])

    if len(alignment_set) == 0:
        return []

    out = run_command(f"LAshow4pathplot -a {db_file} {las_file} {read_id} | sed 's/,//g' | awk 'BEGIN {{first = 1}} NF == 7 {{if (first == 1) {{first = 0}} else {{printf(\"%s\\n%s\\n%s\\n%s\\n\", header, aseq, bseq, symbol)}}; header = $0; aseq = \"\"; bseq = \"\"; symbol = \"\"; count = 0;}} NF < 7 {{if (count == 0) {{aseq = aseq $0}} else if (count == 1) {{symbol = symbol $0}} else {{bseq = bseq $0}}; count++; count %= 3;}} END {{printf(\"%s\\n%s\\n%s\\n%s\\n\", header, aseq, bseq, symbol)}}'").strip().split("\n")

    ret = []   # [(ab, ae, bb, be, FlattenCigar), ...]
    for header, aseq, bseq, symbol in zip(*([iter(out)] * 4)):    # split every 4 lines (= single entry)
        alignment = tuple(map(int, header.replace(" ", "").split("\t")[2:6]))   # ab, ae, bb, be
        if alignment in alignment_set:
            ret.append((*alignment, convert_symbol(*map(lambda x: x[slice(*find_boundary(aseq, bseq))],
                                                        [aseq, bseq, symbol]))))
    logger.debug(f"read {read_id}: {len(ret)} paths loaded")
    return sorted(ret, key=lambda x: x[2])   # sort by bb


def split_tr(ab, bb, fcigar):
    apos, bpos = ab, bb
    unit_intvls = [(bpos, apos)]
    # Iteratively find max{ax} such that bx == (last unit end)
    for i, c in enumerate(fcigar):
        if c != "I":
            apos += 1
        if c != "D":
            bpos += 1
        if bpos == unit_intvls[-1][1] and (i == len(fcigar) - 1 or fcigar[i + 1] != "D"):
            unit_intvls.append((unit_intvls[-1][1], apos))
    return unit_intvls


def find_units_single(read_id, tr_intervals, alignments, db_file, las_file, max_cv=0.1):
    ret = []
    alignment_set = filter_alignments(tr_intervals, alignments)
    for ab, ae, bb, be, fcigar in load_paths(read_id, alignment_set, db_file, las_file):
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


def load_dumps(start_dbid, end_dbid, db_file, las_file):
    # Extract data from DBdump"s and LAdump"s output
    return ({read_id: list(df.apply(lambda d: tuple(d[["start", "end"]]), axis=1))
             for read_id, df
             in (pd.DataFrame([list(map(int, x.split("\t")))
                               for x in run_command(f"DBdump -rh -mtan {db_file} {start_dbid}-{end_dbid} | awk '$1 == \"R\" {{dbid = $2}} $1 == \"T0\" && $2 > 0 {{for (i = 1; i <= $2; i++) printf(\"%s\\t%s\\t%s\\n\", dbid, $(2 * i + 1), $(2 * i + 2))}}'").strip().split("\n")],
                              columns=("dbid", "start", "end")) \
                 .groupby("dbid"))},
            {read_id: list(df.sort_values(by="abpos", kind="mergesort") \
                           .sort_values(by="distance", kind="mergesort") \
                           .apply(lambda d: tuple((*map(int, d[["abpos", "aepos", "bbpos", "bepos", "distance"]]), d["slope"])), axis=1))
             for read_id, df
             in (pd.DataFrame([list(map(int, x.split("\t")))
                               for x in run_command(f"LAdump -c {db_file} {las_file} {start_dbid}-{end_dbid} | awk '$1 == \"P\" {{dbid = $2}} $1 == \"C\" {{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\", dbid, $2, $3, $4, $5)}}'").strip().split("\n")],
                              columns=("dbid", "abpos", "aepos", "bbpos", "bepos")) \
                 .assign(distance=lambda x: x["abpos"] - x["bbpos"]) \
                 .assign(slope=lambda x: ((x["aepos"] - x["abpos"]) / (x["bepos"] - x["bbpos"])).round(3)) \
                 .groupby("dbid"))})


def find_units_multi(start_dbid, end_dbid, db_file, las_file):
    tr_intervals_all, alignments_all = load_dumps(start_dbid, end_dbid, db_file, las_file)
    return [find_units_single(read_id, tr_intervals_all[read_id], alignments_all[read_id], db_file, las_file)
            for read_id in sorted(tr_intervals_all.keys())]


def find_units(start_dbid, end_dbid, db_file, las_file, n_core, out_file):
    """Split the tasks ranging from <start_dbid> and <end_dbid> into <n_core> sub-tasks.
    Each splitted sub-task is then parallely executed, and output into a file.
    A function named "find_units_single" above is actually the core function.
    """
    if n_core == 1:
        rets = find_units_multi(start_dbid, end_dbid, db_file, las_file)
    else:
        unit_n = -(-(end_dbid - start_dbid + 1) // n_core)
        args = [(start_dbid + i * unit_n,   # start_dbid
                 min([start_dbid + (i + 1) * unit_n - 1, end_dbid]),   # end_dbid
                 db_file, las_file)
                for i in range(n_core)]

        with Pool(n_core) as pool:
            rets = [ret for rets in pool.starmap(find_units_multi, args) for ret in rets]

    # Aggregate results into pd.DF
    results = {}
    index = 0
    for ret in rets:
        for r in ret:
            results[index] = r
            index += 1
    pd.DataFrame.from_dict(results, orient="index",
                           columns=("read_id", "start", "end", "mean_ulen", "cv_ulen", "units")) \
                .to_csv(out_file, sep="\t")


def concat_df(dir_name, prefix, sep="\t", index_col=0):
    return pd.concat([pd.read_csv(fname, sep="\t", index_col=index_col)
                      for fname in run_command(f"find {dir_name} -name '{prefix}.*' | sort").strip().split("\n")]) \
             .reset_index(drop=True)


if __name__ == "__main__":
    """Only for internal usage by run_datruf."""
    p = argparse.ArgumentParser()
    p.add_argument("db_fname", type=str)
    p.add_argument("las_fname", type=str)
    p.add_argument("out_fname", type=str)
    p.add_argument("start_dbid", type=int)
    p.add_argument("end_dbid", type=int)
    p.add_argument("n_core", type=int)
    args = p.parse_args()

    find_units(args.start_dbid, args.end_dbid, args.db_fname, args.las_fname, args.n_core, args.out_fname)
