from collections import defaultdict, Counter
from logzero import logger
from BITS.plot.plotly import make_hist, make_layout, show_plot
from BITS.util.union_find import UnionFind


def read_id_to_overlaps(read_id, overlaps):
    return list(filter(lambda o: o.a_read_id == read_id or o.b_read_id == read_id, overlaps))


def reduce_same_overlaps(overlaps, max_bp_slip=100):
    """Reduce overlaps having same positions into single overlap with the best score.
    Overlaps whose difference in bp is at most `max_bp_slip` are regarded as same.
    """
    def are_similar(x, y):
        """Determine whether two overlaps, `x` and `y`, have similar overlapping positions."""
        return (abs(x.a_start - y.a_start) <= max_bp_slip
                and abs(x.a_end - y.a_end) <= max_bp_slip
                and abs(x.b_start - y.b_start) <= max_bp_slip
                and abs(x.b_end - y.b_end) <= max_bp_slip)

    # Aggregate overlaps by read pair
    ovlps_by_pair = defaultdict(list)
    for o in overlaps:
        ovlps_by_pair[(o.a_read_id, o.b_read_id, o.strand)].append(o)
    # For each read pair, aggregate overlaps by position
    reduced_overlaps = []
    for read_pair, ovlps in ovlps_by_pair.items():
        uf = UnionFind(len(ovlps))
        for i in range(len(ovlps)):
            for j in range(len(ovlps)):
                if i >= j:
                    continue
                if uf.in_same_set(i, j):
                    continue
                if are_similar(ovlps[i], ovlps[j]):
                    uf.unite(i, j)
        ovlps_by_pos = {}
        for i in range(len(ovlps)):
            parent = uf.get_root(i)
            if parent not in ovlps_by_pos or ovlps[i].diff < ovlps_by_pos[parent].diff:
                ovlps_by_pos[parent] = ovlps[i]
        reduced_overlaps += list(ovlps_by_pos.values())
    logger.info(f"#Overlaps: {len(overlaps)} -> {len(reduced_overlaps)}")
    return reduced_overlaps


def filter_overlaps(overlaps, max_diff=2., min_ovlp_len=3000):
    """Filter overlaps based on the maximum sequence dissimilarity and minimum overlap length."""
    f = lambda o: ((o.a_end - o.a_start + o.b_end - o.b_start) // 2 >= min_ovlp_len
                   and o.diff < max_diff)
    filtered_overlaps = list(filter(f, overlaps))
    logger.info(f"#Overlaps: {len(overlaps)} -> {len(filtered_overlaps)}")
    return filtered_overlaps


def best_overlaps_per_pair(overlaps):
    """Keep only one overlap for each read pair (+ strand), namely best-overlap logic for
    slippy overlaps."""
    ovlp_by_pair = {}
    for o in overlaps:
        read_pair = (o.a_read_id, o.b_read_id, o.strand)
        if read_pair not in ovlp_by_pair or o.diff < ovlp_by_pair[read_pair].diff:
            ovlp_by_pair[read_pair] = o
    best_overlaps = sorted(ovlp_by_pair.values())
    logger.info(f"#Overlaps: {len(overlaps)} -> {len(best_overlaps)}")
    return best_overlaps


def best_overlaps(overlaps):
    """Best-overlap logic, i.e., keep only one best in-edge and one best out-edge for each read."""
    pass


def plot_n_ovlps_per_read(overlaps, reads):
    unique_overlaps = set([tuple(sorted((o.a_read_id, o.b_read_id))) for o in overlaps])
    n_ovlps = Counter()   # {read_id: ovlp_count}
    for a_read, b_read in unique_overlaps:
        n_ovlps[a_read] += 1
        n_ovlps[b_read] += 1

    stratified_reads = defaultdict(list)   # {total_unit_length_in_kb: reads}
    for read in reads:
        stratified_reads[sum([unit.length for unit in read.units]) // 1000].append(read)

    traces = [make_hist([n_ovlps[read.id] for read in sreads],
                        bin_size=1,
                        name=f"{covered_kb}-{covered_kb + 1}kb-covered reads")
              for covered_kb, sreads in sorted(stratified_reads.items())]
    layout = make_layout(title=(f"Stratified by reads according to the total unit length"),
                             x_title="Number of overlaps", y_title="Frequency")
    layout["barmode"] = "stack"
    show_plot(traces, layout)
