from collections import defaultdict
from matplotlib.colors import XKCD_COLORS
from BITS.plot.plotly import make_line, make_scatter, make_layout, show_plot
from .overlapper.ava_labeled_reads import _ava_labeled_reads


def show_labeled_read(read, max_read_len):
    shapes = [make_line(0, 0, read.length, 0, "grey", 3, 'above')]
    traces = []
    id_to_col = {}
    for i, col in enumerate(XKCD_COLORS.values()):   # NOTE: up to ~1000 cols
        id_to_col[i] = col
    shapes += [make_line(unit.start, 0, unit.end, 0, id_to_col[unit.repr_id], 5, 'above')
               for i, unit in enumerate(read.units)]
    trace_unit = make_scatter([unit.start for i, unit in enumerate(read.units)],
                              [0 for unit in read.units],
                              text=[unit.repr_id for unit in read.units],
                              text_pos="top right", text_size=10, text_col="black", mode="text",
                              name="TR unit",
                              show_legend=False)
    trace_info = make_scatter([unit.start for i, unit in enumerate(read.units)],
                              [0 for unit in read.units],
                              text=[f"unit {i} (id={unit.repr_id if read.synchronized else ''})<br>"
                                    f"[{unit.start}:{unit.end}] ({unit.length} bp)"
                                    for i, unit in enumerate(read.units)],
                              col=[id_to_col[unit.repr_id] for unit in read.units],
                              show_legend=False)
    traces += [trace_unit, trace_info]

    # Show plot
    layout = make_layout(height=200, title=f"{read.id}",
                         x_range=(-0.01 * max_read_len, max_read_len * 1.01),
                         x_grid=False, y_grid=False, shapes=shapes)
    show_plot(traces, layout)


def show_labeled_reads(labeled_reads):
    max_length = max([read.length for read in labeled_reads])
    for read in labeled_reads:
        show_labeled_read(read, max_length)


def show_labeled_reads_pileup(read_id, labeled_reads_by_id, height=None):
    labeled_reads = labeled_reads_by_id[read_id]
    labeled_reads_by_id_single = {read.id: read for read in labeled_reads}

    rel_start_pos = defaultdict(list)
    rel_start_pos[read_id].append(0)

    for o in _ava_labeled_reads(labeled_reads):
        if o.a_read_id == read_id:
            rel_start_pos[o.b_read_id].append(o.a_start - o.b_start)
        elif o.b_read_id == read_id:
            rel_start_pos[o.a_read_id].append(o.b_start - o.a_start if o.strand == 0
                                              else (o.b_len - o.b_end) - (o.a_len - o.a_end))

    x_min = min([min(poss) for read_id, poss in rel_start_pos.items()])
    x_max = max([max(poss) + labeled_reads_by_id_single[read_id].length
                 for read_id, poss in rel_start_pos.items()])

    start_pos_to_read = sorted([(pos, read_id)
                                for read_id, poss in rel_start_pos.items()
                                for pos in poss])

    y = 0
    traces, shapes = [], []
    for start_pos, read_id in start_pos_to_read:
        read = labeled_reads_by_id_single[read_id]
        shapes += [make_line(start_pos, y, start_pos + read.length, y, "grey", 3, 'above')]
        id_to_col = {}
        for i, col in enumerate(XKCD_COLORS.values()):   # NOTE: up to ~1000 cols
            id_to_col[i] = col
        shapes += [make_line(start_pos + unit.start, y, start_pos + unit.end, y,
                             id_to_col[unit.repr_id], 5, 'above')
                   for i, unit in enumerate(read.units)]
        trace_unit = make_scatter([start_pos + unit.start for i, unit in enumerate(read.units)],
                                  [y for unit in read.units],
                                  show_legend=False)
        trace_info = make_scatter([start_pos + unit.start for i, unit in enumerate(read.units)],
                                  [y for unit in read.units],
                                  text=[f"unit {i} (id={unit.repr_id if read.synchronized else ''})<br>"
                                        f"[{unit.start}:{unit.end}] ({unit.length} bp)"
                                        for i, unit in enumerate(read.units)],
                                  col=[id_to_col[unit.repr_id] for unit in read.units],
                                  show_legend=False)
        traces += [trace_unit, trace_info]
        y -= 1

    traces += [make_scatter([start_pos for start_pos, read_id in start_pos_to_read],
                            [-y for y in range(len(start_pos_to_read))],
                            text=[read_id for start_pos, read_id in start_pos_to_read],
                            col="grey", marker_size=1, show_legend=False)]

    layout = make_layout(height=100 + 50 * len(start_pos_to_read) if height is None else height,
                         x_range=(x_min, x_max), x_grid=False, y_grid=False, x_zeroline=False,
                         shapes=shapes)
    show_plot(traces, layout)
