from BITS.seq.align import EdlibRunner
from BITS.seq.utils import reverse_seq

er_prefix = EdlibRunner("prefix", revcomp=False, cyclic=False)


def can_be_query(focal_seq, opponent_seq):
    """Check if `focal_seq` can be a query for `opponent_seq`. If `focal_seq` is too long to
    map to `opponent_seq`, return False.
    """
    if len(focal_seq) < len(opponent_seq) * 1.1:
        return True
    return False


def prefix_alignment(query, target):
    """Compute prefix alignment between `query` and `target`. That is, start positions of the
    alignment are 0 for both sequences, but end positions are not constrained.
    """
    assert len(query) > 0 and len(target) > 0, "Empty sequence is not allowed"
    aln = None
    if can_be_query(query, target):
        aln = er_prefix.align(query, target)   # map `query` to `target`
        q_end, t_end = len(query), aln.t_end
    if can_be_query(target, query):
        aln_swap = er_prefix.align(target, query)   # map `target` to `query`
        if aln is None or aln.diff > aln_swap.diff:
            aln = aln_swap
            q_end, t_end = aln.t_end, len(target)
    assert aln is not None, "Both sequences were not query"
    return (aln, q_end, t_end)


def suffix_alignment(query, target):
    return prefix_alignment(reverse_seq(query), reverse_seq(target))


def dovetail_alignment(query, target, q_match_pos, t_match_pos):
    """Compute dovetail alignment between `query` and `target` given positions which
    confidently match between them. The alignment will be splitted into the following parts:
      1. Best suffix alignment between `query[:q_match_pos]` and `target[:t_match_pos]`
      2. Best prefix alignment between `query[q_match_pos:]` and `target[t_match_pos:]`
    """
    assert 0 <= q_match_pos <= len(query), f"`q_match_pos` out of range"
    assert 0 <= t_match_pos <= len(target), f"`t_match_pos` out of range"
    
    aln_len_tot, aln_n_diff_tot = 0, 0

    # Alignment up to `[q|t]_match_pos`
    if q_match_pos == 0 or t_match_pos == 0:
        q_start, t_start = q_match_pos, t_match_pos
    else:
        aln_first, q_first, t_first = suffix_alignment(query[:q_match_pos], target[:t_match_pos])
        q_start, t_start = q_match_pos - q_first, t_match_pos - t_first
        aln_len_tot += aln_first.length
        aln_n_diff_tot += int(aln_first.length * aln_first.diff)

    # Alignment from `[q|t]_match_pos`
    if q_match_pos == len(query) or t_match_pos == len(target):
        q_end, t_end = q_match_pos, t_match_pos
    else:
        aln_second, q_second, t_second = prefix_alignment(query[q_match_pos:], target[t_match_pos:])
        q_end, t_end = q_match_pos + q_second, t_match_pos + t_second
        aln_len_tot += aln_second.length
        aln_n_diff_tot += int(aln_second.length * aln_second.diff)

    return (q_start, q_end, t_start, t_end, aln_len_tot, aln_n_diff_tot / aln_len_tot)
