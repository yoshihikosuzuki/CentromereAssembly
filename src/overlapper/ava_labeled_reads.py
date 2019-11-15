from multiprocessing import Pool
from logzero import logger
from BITS.seq.align import EdlibRunner
from vca.overlapper.dovetail_overlap import dovetail_alignment
from vca.types import Overlap

er_global = EdlibRunner("global", revcomp=False, cyclic=False)


def svs_labeled_reads(a_read, b_read, repr_alignments,
                      min_n_units=3, max_units_diff=0.01, max_seq_diff=0.02):
    """Overlap by sequence identity only on representative units. Sequence dovetail overlap including
    non-TR regions is done just for excluding overlaps with too different non-TR regions."""
    overlaps = set()
    if min(len(a_read.units), len(b_read.units)) < min_n_units:
        return overlaps

    # from  a ----->     to  a ----->
    #       b    ----->      b ----->
    b_start_unit = 0
    for a_start_unit in range(len(a_read.units) - min_n_units + 1):
        alignments = [repr_alignments[(a_read.units[a_start_unit + i].repr_id,
                                       b_read.units[b_start_unit + i].repr_id)]
                      for i in range(min(len(a_read.units) - a_start_unit, len(b_read.units)))]
        diff = sum([a.length * a.diff for a in alignments]) / sum([a.length for a in alignments])
        if diff < max_units_diff:
            # Confirm sequences including non-TR regions are not so much different
            overlap = dovetail_alignment(a_read.seq, b_read.seq,
                                         a_read.units[a_start_unit].start,
                                         b_read.units[b_start_unit].start)
            a_start, a_end, b_start, b_end, seq_len, seq_diff = overlap
            if seq_diff < max_seq_diff:
                if a_read.strand == 1:
                    a_start, a_end = a_read.length - a_end, a_read.length - a_start
                    b_start, b_end = b_read.length - b_end, b_read.length - b_start
                overlaps.add(Overlap(a_read.id, b_read.id, 0 if a_read.strand == b_read.strand else 1,
                                     a_start, a_end, a_read.length, b_start, b_end, b_read.length,
                                     round(100 * diff, 2)))
    
    # from  a  ----->  to     a ----->   
    #       b ----->      b ----->
    a_start_unit = 0
    for b_start_unit in range(1, len(b_read.units) - min_n_units + 1):
        alignments = [repr_alignments[(a_read.units[a_start_unit + i].repr_id,
                                       b_read.units[b_start_unit + i].repr_id)]
                      for i in range(min(len(b_read.units) - b_start_unit, len(a_read.units)))]
        diff = sum([a.length * a.diff for a in alignments]) / sum([a.length for a in alignments])
        if diff < max_units_diff:
            # Confirm sequences including non-TR regions are not so much different
            overlap = dovetail_alignment(a_read.seq, b_read.seq,
                                         a_read.units[a_start_unit].start,
                                         b_read.units[b_start_unit].start)
            a_start, a_end, b_start, b_end, seq_len, seq_diff = overlap
            if seq_diff < max_seq_diff:
                if a_read.strand == 1:
                    a_start, a_end = a_read.length - a_end, a_read.length - a_start
                    b_start, b_end = b_read.length - b_end, b_read.length - b_start
                overlaps.add(Overlap(a_read.id, b_read.id, 0 if a_read.strand == b_read.strand else 1,
                                     a_start, a_end, a_read.length, b_start, b_end, b_read.length,
                                     round(100 * diff, 2)))
    
    return overlaps


def _ava_labeled_reads(labeled_reads):
    # Precompute all-vs-all global alignments between the representative units
    for i in range(len(labeled_reads) - 1):
        assert labeled_reads[i].repr_units == labeled_reads[i + 1].repr_units, \
            "Representative units of the reads must be same"
    repr_units = labeled_reads[0].repr_units
    repr_alignments = {}
    for id_i, seq_i in repr_units.items():
        for id_j, seq_j in repr_units.items():
            if id_i > id_j:
                continue
            repr_alignments[(id_i, id_j)] = repr_alignments[(id_j, id_i)] = er_global.align(seq_i, seq_j)
                
    # Compute all-vs-all labeled_read overlap
    overlaps = set()
    for read_i in labeled_reads:
        for read_j in labeled_reads:
            if read_i.id >= read_j.id:
                continue
            overlaps.update(svs_labeled_reads(read_i, read_j, repr_alignments))
    return overlaps


def ava_labeled_reads(labeled_reads_by_id, n_core=1):
    overlaps = set()
    with Pool(n_core) as pool:
        for ovlps in pool.map(_ava_labeled_reads, list(labeled_reads_by_id.values())):
            overlaps.update(ovlps)
    return sorted(overlaps)
