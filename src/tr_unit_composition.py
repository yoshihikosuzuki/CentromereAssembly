from collections import Counter
from BITS.plot.plotly import make_scatter, make_layout, show_plot


def read_to_ulen_comp(read):
    return Counter([unit.length for unit in read.units]).most_common()   # [(ulen, count)]


def scatter_ulen_comp(read):
    """Draw a unit length vs unit length * count plot for a single read."""
    comp = sorted(read_to_ulen_comp(read))
    show_plot([make_scatter([ulen for ulen, count in comp],
                            [ulen * count for ulen, count in comp],
                            mode="lines+markers", marker_size=8, show_legend=False)],
              make_layout(x_title="Unit length [bp]",
                          y_title="Total length of the units of the length [bp]"))


def scatter_ulen_comps(reads):
    """Draw a unit length vs unit length * count plot for multiple reads."""
    comps = [sorted(read_to_ulen_comp(read)) for read in reads]
    show_plot([make_scatter([ulen for ulen, count in comp],
                            [ulen * count for ulen, count in comp],
                            mode="lines+markers", marker_size=6, name=f"read {reads[i].id}")
               for i, comp in enumerate(comps)],
              make_layout(x_title="Unit length [bp]",
                          y_title="Total length of the units of the length [bp]"))


def show_ulen_dist(read):
    """Draw a positional distribution of unit lengths on `read`."""
    show_plot([make_scatter([unit.start for unit in read.units],
                            [unit.length for unit in read.units],
                            show_legend=False)],
               make_layout(x_title="Start position",
                           y_title="Unit length [bp]"))
