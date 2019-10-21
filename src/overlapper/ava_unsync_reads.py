import argparse
from dataclasses import dataclass
import numpy as np
from logzero import logger
from BITS.seq.align import EdlibRunner
from BITS.seq.utils import reverse_seq
from BITS.util.io import save_pickle, load_pickle
from BITS.util.proc import NoDaemonPool, run_command
from BITS.util.scheduler import Scheduler
from vca.types import Overlap, revcomp_read

out_dir        = "ava_unsync"
out_prefix     = "ovlps"
scatter_prefix = "run_ava_unsync"
gather_fname   = f"{out_dir}/gather.sh"
log_fname      = f"{out_dir}/log"


@dataclass(eq=False)
class UnsyncReadsOverlapper:
    """Class for executing all-vs-all overlap between unsynchronized TR reads.

    positional arguments:
      @ n_distribute <int> : Number of jobs
      @ n_core       <int> : Number of cores per job

    optional arguments:
      @ scheduler              <Scheduler> [Scheduler("sge", "qsub", "all.q")]
          : Job scheduler
      @ centromere_reads_fname <str>       ["centromere_reads.pkl"]
          : File of centromere reads
      @ out_fname              <str>       ["centromere_reads_unsync_overlaps.pkl"]
          : Output file name
    """
    n_distribute          : int
    n_core                : int
    scheduler             : Scheduler = Scheduler("sge", "qsub", "all.q")
    centromere_reads_fname: str       = "centromere_reads.pkl"
    out_fname             : str       = "centromere_reads_unsync_overlaps.pkl"

    def __post_init__(self):
        run_command(f"mkdir -p {out_dir}; rm -f {out_dir}/*")

    def run(self):
        jids = []
        for i in range(self.n_distribute):
            index = str(i + 1).zfill(int(np.log10(self.n_distribute) + 1))
            out_fname = f"{out_dir}/{out_prefix}.{index}.pkl"
            script_fname = f"{out_dir}/{scatter_prefix}.{index}.sh"
            script = (f"python -m vca.ava_unsync_reads {self.centromere_reads_fname} "
                      f"{out_fname} {self.n_distribute} {self.n_core} {i}")
            
            jids.append(self.scheduler.submit(script,
                                              script_fname,
                                              job_name="ava_unsync",
                                              log_fname=log_fname,
                                              n_core=self.n_core))

        self.scheduler.submit("sleep 1s",
                              gather_fname,
                              job_name="ava_unsync_gather",
                              log_fname=log_fname,
                              depend=jids,
                              wait=True)

        merged = []
        fnames = run_command(f"find {out_dir} -name '{out_prefix}.*' | sort").strip().split('\n')
        for fname in fnames:
            merged += load_pickle(fname)
        save_pickle(sorted(merged), self.out_fname)


def seq_to_spectrum(seq, k=13):
    return set([seq[i:i + k] for i in range(len(seq) - k + 1)])


def read_to_boundary_k_monomers(read, k=2, offset=1):
    """Extract <k>-monomers (including sequences between monomers) of the boundaries.
    Sequences of both orientations are returned. Therefore, the total number of <k>-monomers are 4.
    <offset> boundary monomers are skipped to avoid noise on sequence."""
    if len(read.units) < 2 * (k + offset):
        logger.info(f"Read {read.id}: # of units not sufficient. Skip.")
        return []

    prefix_k_unit = read.seq[read.units[offset].start:read.units[offset + k - 1].end]
    suffix_k_unit = read.seq[read.units[-offset - k].start:read.units[-offset - 1].end]
    return ((prefix_k_unit, offset, seq_to_spectrum(prefix_k_unit)),
            (suffix_k_unit, len(read.units) - offset - k, seq_to_spectrum(suffix_k_unit)))


def reads_to_boundary_k_monomers(reads, k=2, offset=1):
    return {(read.id, strand): read_to_boundary_k_monomers(read if strand == 0 else revcomp_read(read))
            for read in centromere_reads for strand in (0, 1)}


er_global = EdlibRunner("global", revcomp=False, cyclic=False)
er_glocal = EdlibRunner("glocal", revcomp=False, cyclic=False)
er_prefix = EdlibRunner("local", revcomp=False, cyclic=False)


def _calc_dovetail_alignment(a_read, b_read, a_index, b_index, k, max_init_diff=0.02):
    """Compute dovetail alignment between two reads starting by mapping `k`-units of `b_read` to (`k+1`)-units of `a_read`.
    If the initla unit global alignment
    """
    a_seq = a_read.seq[a_read.units[a_index].start:a_read.units[a_index + k].end]
    b_seq = b_read.seq[b_read.units[b_index].start:b_read.units[b_index + k - 1].end]

    alignment = er_glocal.align(b_seq, a_seq)
    if alignment.diff > max_init_diff:   # initial screening
        return (-1, -1, 100.)
    a_start = a_read.units[a_index].start + alignment.t_start
    b_start = b_read.units[b_index].start

    tot_alignment_length = 0
    tot_alignment_diff = 0

    # Compute forward dovetail overlap staring at (ab, bb) = (a_start, b_start)
    logger.debug("###")
    logger.debug(f"a_units[{a_index}-{a_index + k}] ({a_read.id}[{a_read.units[a_index].start}:{a_read.units[a_index + k].end}]) vs "
                 f"b_unit[{b_index}-{b_index + k - 1}] ({b_read.id}[{b_read.units[b_index].start}:{b_read.units[b_index + k - 1].end}])")
    logger.debug(f"a_start @ {a_start} ~ b_start @ {b_start} "
                 f"({round(100 * alignment.diff, 1)} %diff, {alignment.length} bp)")

    # First, up to the start position
    if a_start == 0 and b_start == 0:
        a_start_pos, b_start_pos = 0, 0
    else:
        a_seq = reverse_seq(a_read.seq[:a_start])   # reverse sequences so that prefix alignment can be taken
        b_seq = reverse_seq(b_read.seq[:b_start])
        query, target = f"rev(a[0:{a_start}])", f"rev(b[0:{b_start}])"
        prefix_alignment = None
        if len(a_seq) < len(b_seq) * (1 + max_init_diff):   # a_seq can be query
            prefix_alignment = er_prefix.align(a_seq, b_seq)
            a_start_pos, b_start_pos = 0, b_start - prefix_alignment.t_end
        if prefix_alignment is None or len(b_seq) < len(a_seq) * (1 + max_init_diff):   # b_seq can be query; both can hold simultaneously
            swap_prefix_alignment = er_prefix.align(b_seq, a_seq)
            if prefix_alignment is None or prefix_alignment.diff > swap_prefix_alignment.diff:
                query, target = target, query
                prefix_alignment = swap_prefix_alignment
                a_start_pos, b_start_pos = a_start - prefix_alignment.t_end, 0
        tot_alignment_length += prefix_alignment.length
        tot_alignment_diff += int(prefix_alignment.diff * prefix_alignment.length)
        logger.debug(f"prefix {query} -> {target}[{prefix_alignment.t_start}:{prefix_alignment.t_end}] "
                     f"({round(100 * prefix_alignment.diff, 1)} %diff, {prefix_alignment.length} bp)")

    # Second, after the start position
    a_seq = a_read.seq[a_start:]
    b_seq = b_read.seq[b_start:]
    query, target = f"a[{a_start}:{a_read.length}]", f"b[{b_start}:{b_read.length}]"
    suffix_alignment = None
    if len(a_seq) < len(b_seq) * (1 + max_init_diff):   # a_seq can be query
        suffix_alignment = er_prefix.align(a_seq, b_seq)
        a_end_pos, b_end_pos = a_read.length, b_start + suffix_alignment.t_end
    if suffix_alignment is None or len(b_seq) < len(a_seq) * (1 + max_init_diff):   # b_seq can be query; both can hold simultaneously
        swap_suffix_alignment = er_prefix.align(b_seq, a_seq)
        if suffix_alignment is None or suffix_alignment.diff > swap_suffix_alignment.diff:
            query, target = target, query
            suffix_alignment = swap_suffix_alignment
            a_end_pos, b_end_pos = a_start + suffix_alignment.t_end, b_read.length
    tot_alignment_length += suffix_alignment.length
    tot_alignment_diff += int(suffix_alignment.diff * suffix_alignment.length)
    logger.debug(f"suffix {query} -> {target}[{suffix_alignment.t_start}:{suffix_alignment.t_end}] "
                 f"({round(100 * suffix_alignment.diff, 1)} %diff, {suffix_alignment.length} bp)")

    tot_diff = tot_alignment_diff / tot_alignment_length
    logger.debug(f"{round(100 * tot_diff, 1)} %diff, {tot_alignment_length} bp")
    return (a_start_pos, a_end_pos, b_start_pos, b_end_pos, tot_diff)


def align_b_to_a(a_read, b_read, b_index, k):
    return [_calc_dovetail_alignment(a_read, b_read, a_index, b_index, k)   # k+1-unit for a; k-unit for b
            for a_index in range(len(a_read.units) - k)]


def sva_overlap(a_read, centromere_reads_by_id, read_specs, boundary_k_units, k=2, max_diff=0.02):
    overlaps = set()
    if len(a_read.units) < k + 1:
        return overlaps

    for (b_read_id, strand), k_units in boundary_k_units.items():
        if a_read.id == b_read_id:
            continue
        if len(k_units) == 0:
            continue
        for k_unit_seq, k_unit_start, k_unit_spec in k_units:   # TODO: k_unit_seq も不要
            if len(read_specs[a_read.id] & k_unit_spec) / len(k_unit_spec) < 0.1:
                continue
            b_read = centromere_reads_by_id[b_read_id]
            if strand == 1:
                b_read = revcomp_read(b_read)
            # map b_read.seq[b_read.units[k_unit_start].start:b_read.units[k_unit_start + k - 1].end] to a_read.seq
            for a_start, a_end, b_start, b_end, diff in align_b_to_a(a_read, b_read, k_unit_start, k):
                if diff < max_diff:
                    overlaps.add(Overlap(a_read.id, b_read_id, strand,
                                         a_start, a_end, a_read.length,
                                         b_start, b_end, b_read.length, diff))
    return overlaps


if __name__ == "__main__":
    """Only for internal usage."""
    p = argparse.ArgumentParser()
    p.add_argument("centromere_reads_fname", type=str)
    p.add_argument("out_fname", type=str)
    p.add_argument("n_distribute", type=int)
    p.add_argument("n_core", type=int)
    p.add_argument("index", type=int)
    args = p.parse_args()

    centromere_reads = load_pickle(args.centromere_reads_fname)
    centromere_reads_by_id = {read.id: read for read in centromere_reads}
    # K-mer spectrum of the reads
    read_specs = {read.id: seq_to_spectrum(read.seq) for read in centromere_reads}
    # Boundary units from all the reads
    boundary_k_units = reads_to_boundary_k_monomers(centromere_reads)

    overlaps = set()
    unit_n = -(-len(centromere_reads) // args.n_distribute)
    with NoDaemonPool(args.n_core) as pool:
        for ret in pool.starmap(sva_overlap,
                                [(read, centromere_reads_by_id, read_specs, boundary_k_units)
                                 for read in centromere_reads[args.index * unit_n:(args.index + 1) * unit_n]]):
            overlaps.update(ret)

    save_pickle(overlaps, args.out_fname)
