from logzero import logger
from BITS.seq.align import EdlibRunner
from vca.overlapper.dovetail_overlap import dovetail_alignment
from vca.types import Overlap

er_glocal = EdlibRunner("glocal", revcomp=False, cyclic=False)


def _svs_overlap_forward(boundary_read, whole_read, boundary_start, boundary_end,
                         k_for_unit, min_kmer_ovlp, max_init_diff,
                         read_forward_specs, read_boundary_specs):
    match_poss = set()
    
    # Filter by k-mer spectrum
    boundary_spec = read_boundary_specs[(boundary_read.id, boundary_read.strand,
                                         boundary_start, boundary_end)]
    whole_spec = read_forward_specs[whole_read.id]
    if len(whole_spec & boundary_spec) / len(boundary_spec) < min_kmer_ovlp:
        return match_poss

    # Filter by mapping of k-unit
    boundary_seq = boundary_read.seq[boundary_start:boundary_end]
    aln = er_glocal.align(boundary_seq, whole_read.seq)
    if aln.diff > max_init_diff:
        return match_poss
    logger.debug(f"{boundary_read.id}{'' if boundary_read.strand == 0 else '*'}"
                 f"[{boundary_start}:{boundary_end}]"
                 f" -> {whole_read.id}[{aln.t_start}:{aln.t_end}]")

    # For each (k+1)-unit of `whole_read`, map k-unit of `boundary_read`
    for i in range(len(whole_read.units) - k_for_unit):
        whole_seq = whole_read.seq[whole_read.units[i].start:whole_read.units[i + k_for_unit].end]
        aln = er_glocal.align(boundary_seq, whole_seq)
        if aln.diff > max_init_diff:
            continue
        boundary_match_pos = boundary_start
        whole_match_pos = whole_read.units[i].start + aln.t_start
        match_poss.add((boundary_match_pos, whole_match_pos))

    if len(match_poss) == 0:
        logger.warning(f"k-unit of read {boundary_read.id} is mapped to "
                       f"whole read {whole_read.id} but not to (k+1)-unit")
    return match_poss


def svs_overlap_forward(boundary_read, whole_read, offset, k_for_unit, min_kmer_ovlp, max_init_diff,
                        read_forward_specs, read_boundary_specs):
    assert whole_read.strand == 0, "`whole_read` must be forward"
    # prefix boundary k-units of `boundary_read` vs `whole_read`
    match_poss = _svs_overlap_forward(boundary_read, whole_read,
                                      boundary_read.units[offset].start,
                                      boundary_read.units[offset + k_for_unit - 1].end,
                                      k_for_unit, min_kmer_ovlp, max_init_diff)
    # suffix boundary
    match_poss.update(_svs_overlap_forward(boundary_read, whole_read,
                                           boundary_read.units[-offset - k_for_unit].start,
                                           boundary_read.units[-offset - 1].end,
                                           k_for_unit, min_kmer_ovlp, max_init_diff))
    return match_poss


def svs_overlap(a_read, b_read, a_read_rc, b_read_rc,
                offset, k_for_unit, min_kmer_ovlp, max_init_diff, max_diff,
                read_forward_specs, read_boundary_specs):
    match_pos_a_to_b = set([(a_match_pos, b_match_pos, 0)
                            for a_match_pos, b_match_pos in svs_overlap_forward(a_read, b_read)])
    match_pos_b_to_a = set([(a_match_pos, b_match_pos, 0)
                            for b_match_pos, a_match_pos in svs_overlap_forward(b_read, a_read)])
    match_pos_ar_to_b = set([(a_read.length - a_match_pos, b_read.length - b_match_pos, 1)
                             for a_match_pos, b_match_pos in svs_overlap_forward(a_read_rc, b_read)])
    match_pos_br_to_a = set([(a_match_pos, b_match_pos, 1)
                             for b_match_pos, a_match_pos in svs_overlap_forward(b_read_rc, a_read)])
    match_poss = match_pos_a_to_b | match_pos_b_to_a | match_pos_ar_to_b | match_pos_br_to_a

    overlaps = set()
    for a_match_pos, b_match_pos, strand in match_poss:
        a_start, a_end, b_start, b_end, length, diff = \
            dovetail_alignment(a_read.seq, (b_read if strand == 0 else b_read_rc).seq,
                               a_match_pos, b_match_pos)
        if diff > max_diff:
            continue
        overlaps.add(Overlap(a_read.id, b_read.id, strand,
                             a_start, a_end, a_read.length,
                             b_start, b_end, b_read.length,
                             round(100 * diff, 2)))
    return sorted(overlaps)
