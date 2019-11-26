from collections import Counter
from BITS.plot.plotly import make_hist, make_scatter, make_layout, show_plot


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


def plot_ulen_transition(read, out_fname=None):
    show_plot([make_scatter([unit.start for unit in read.units],
                            [unit.length for unit in read.units],
                            show_legend=False)],
               make_layout(title=f"Read {read.id} (strand={read.strand})",
                           x_title="Start position on the read",
                           y_title="Unit length [bp]",
                           x_range=(0, read.length)),
              out_fname=out_fname)


def reads_to_ulens_count(reads):
    return Counter([unit.length for read in reads for unit in read.units])


def reads_to_ulens_tot(reads):
    ulens = reads_to_ulens_count(reads)
    return sorted([(ulen, ulen * count) for ulen, count in ulens.items()])


def plot_ulens_count(reads, min_ulen=50):
    show_plot(make_hist([unit.length for read in reads for unit in read.units
                         if unit.length >= min_ulen],
                        bin_size=1, show_legend=False),
              make_layout(title=f"Unit count for each unit length (>50 bp unit)",
                          x_title="Unit length [bp]",
                          y_title="Unit count",
                          x_range=(1, None)))


def plot_ulens_tot(reads, x_range=(1, None), out_fname=None):
    ulens = reads_to_ulens_tot(reads)
    show_plot(make_scatter([x[0] for x in ulens],
                           [x[1] for x in ulens],
                           mode="lines", show_legend=False),
              make_layout(title="Total length for each unit length",
                          x_title="Unit length [bp]",
                          y_title="Unit length * unit count [bp]",
                          x_range=x_range), out_fname=out_fname)
