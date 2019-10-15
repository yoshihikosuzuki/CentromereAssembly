from logzero import logger
from BITS.seq.align import EdlibRunner
from BITS.seq.utils import revcomp
from BITS.util.io import save_pickle, load_pickle
from BITS.util.proc import NoDaemonPool
from vca.types import TRUnit, TRRead


def seq_to_spectrum(seq, k=13):
    return set([seq[i:i + k] for i in range(len(seq) - k + 1)])


def revcomp_read(read):
    """Return reverse complement of <read> as a new object. <trs> and <alignments> are not copied."""
    return TRRead(seq=revcomp(read.seq), id=read.id, name=read.name,
                  units=[TRUnit(start=read.length - unit.end,
                                end=read.length - unit.start,
                                repr_id=unit.repr_id,
                                strand=(None if unit.strand is None else 1 - unit.strand))
                         for unit in reversed(read.units)],
                  synchronized=read.synchronized, repr_units=read.repr_units)


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


def reverse_seq(seq):
    return ''.join(list(reversed(seq)))


def _calc_dovetail_alignment(a_read, b_read, a_index, b_index, k, max_init_diff=0.02):
    a_seq = a_read.seq[a_read.units[a_index].start:a_read.units[a_index + k].end]
    b_seq = b_read.seq[b_read.units[b_index].start:b_read.units[b_index + k - 1].end]

    alignment = er_glocal.align(b_seq, a_seq)
    if alignment.diff > max_init_diff:
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

    if a_start > 0 and b_start > 0:
        # First, up to the start position
        a_seq = reverse_seq(a_read.seq[:a_start])   # reverse sequences so that prefix alignment can be taken
        b_seq = reverse_seq(b_read.seq[:b_start])
        query, target = f"rev(a[0:{a_start}])", f"rev(b[0:{b_start}])"
        prefix_alignment = None
        if len(a_seq) < len(b_seq) * (1 + max_init_diff):   # a_seq can be query
            prefix_alignment = er_prefix.align(a_seq, b_seq)
        if prefix_alignment is None or len(b_seq) < len(a_seq) * (1 + max_init_diff):   # b_seq can be query; both can hold simultaneously
            swap_prefix_alignment = er_prefix.align(b_seq, a_seq)
            if prefix_alignment is None or prefix_alignment.diff > swap_prefix_alignment.diff:
                prefix_alignment = swap_prefix_alignment
                query, target = target, query
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
    if suffix_alignment is None or len(b_seq) < len(a_seq) * (1 + max_init_diff):   # b_seq can be query; both can hold simultaneously
        swap_suffix_alignment = er_prefix.align(b_seq, a_seq)
        if suffix_alignment is None or suffix_alignment.diff > swap_suffix_alignment.diff:
            suffix_alignment = swap_suffix_alignment
            query, target = target, query
    tot_alignment_length += suffix_alignment.length
    tot_alignment_diff += int(suffix_alignment.diff * suffix_alignment.length)
    logger.debug(f"suffix {query} -> {target}[{suffix_alignment.t_start}:{suffix_alignment.t_end}] "
                 f"({round(100 * suffix_alignment.diff, 1)} %diff, {suffix_alignment.length} bp)")

    tot_diff = tot_alignment_diff / tot_alignment_length
    logger.debug(f"{round(100 * tot_diff, 1)} %diff, {tot_alignment_length} bp")
    return (a_start, b_start, tot_diff)


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
            for a_start, b_start, diff in align_b_to_a(a_read, b_read, k_unit_start, k):
                if diff < max_diff:
                    overlaps.add((a_read.id, b_read_id, strand, a_start, b_start, diff))
    return overlaps


if __name__ == "__main__":
    import sys
    i, n_distribute, n_core, out_fname = sys.argv[1:]
    i, n_distribute, n_core = list(map(int, [i, n_distribute, n_core]))

    centromere_reads = load_pickle("centromere_reads.pkl")
    centromere_reads_by_id = {read.id: read for read in centromere_reads}
    read_specs = {read.id: seq_to_spectrum(read.seq) for read in centromere_reads}

    # Boundary units from all the reads
    boundary_k_units = reads_to_boundary_k_monomers(centromere_reads)

    overlaps = set()
    unit_n = -(-len(centromere_reads) // n_distribute)
    with NoDaemonPool(n_core) as pool:
        for ret in pool.starmap(sva_overlap,
                                [(read, centromere_reads_by_id, read_specs, boundary_k_units)
                                 for read in centromere_reads[i * unit_n:(i + 1) * unit_n]]):
            overlaps.update(ret)

    save_pickle(overlaps, out_fname)
