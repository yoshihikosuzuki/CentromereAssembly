from dataclasses import dataclass, field
from typing import List
from collections import defaultdict
import numpy as np
from logzero import logger
from interval import interval
from sklearn.neighbors import KernelDensity
from BITS.plot.plotly import make_line, make_rect, make_hist, make_scatter, make_layout, show_plot


def read_to_ulens(read, min_ulen, max_ulen):
    """Extract units of [`min_ulen`..`max_ulen`] bp length from `read.units`."""
    return [unit.length for unit in read.units
            if min_ulen <= unit.length <= max_ulen]


def read_to_ulens_in_intvls(read, intvls):
    """Extract units within `intvls` from `read.units`."""
    return [unit.length for unit in read.units
            if unit.length in intvls]


def stratify_by_covered_length(reads, min_ulen, max_ulen):
    """Split `reads` into `List[reads]` according to the total unit length for every 1000 bp."""
    stratified_reads = defaultdict(list)   # {total_unit_length_in_kb: reads}
    for read in reads:
        stratified_reads[sum(read_to_ulens(read, min_ulen, max_ulen)) // 1000].append(read)
    return stratified_reads


@dataclass(eq=False)
class TRReadFilter:
    """Class for extracting putative centromeric reads from `List[TRRead]` by using an assumption
    that centromeric TR units must be abundant.

    Before executing `run()`, It is recomended to find the best parameters by looking at unit length
    distribution with `hist_unit_lengths()` via Jupyter Notebook.

    Workflow example:
      > f = TRReadFilter(your_params)
      > f.hist_unit_lengths(tr_reads, x_min, x_max)   # `tr_reads` is output of datruf
      > # Modify `[min|max]_ulen` and `min_covered_length` here as you like
      > f.find_peak_ulens(tr_reads)
      > # Modify `band_width`, `min_density`, `deviation` here as you like
      > centromere_reads = f.run(tr_reads, show_density=False)

    optional arguments:
      @ min_ulen           <int>   [50]
          : Peaks of the unit length are detected from the range [`min_ulen`..`max_ulen`].
      @ max_ulen           <int>   [500]
      @ min_covered_length <int>   [2000]
          : Reads not covered more than `min_covered_length` bp by `min_ulen`-`max_ulen` bp units
          : are not extracted.
      @ band_width         <int>   [5]
          : Parameter for peak detection. Smaller value results in more sensitive to peaks.
      @ min_density        <float> [0.005]
          : Minimum density (= relative frequency) required for peak unit lengths.
      @ deviation          <float> [0.1]
          : Units inside the range [`peak_ulen * (1 - deviation)`..`peak_ulen * (1 + deviation)`]
            are extracted.
    """
    min_ulen           : int       = 300
    max_ulen           : int       = 400
    min_covered_length : int       = 2000
    band_width         : float     = 0.5
    min_density        : float     = 0.01
    deviation          : float     = 0.015
    peak_ulens         : List[int] = field(init=False, default=None)
    peak_intvls        : interval  = field(init=False, default=None)

    def __post_init__(self):
        if self.min_ulen < 30:
            logger.warn(f"The value `min_ulen={self.min_ulen}` is very small and "
                        f"the calculation might be stacked.")

    def run(self, tr_reads):
        # Find peak unit length(s) and interval(s)
        self.find_peak_ulens(tr_reads, show_density=False)

        # Extract reads using `peak_intvls`
        centromere_reads = \
            list(filter(lambda read: \
                        sum(read_to_ulens_in_intvls(read, self.peak_intvls)) >= self.min_covered_length,
                        tr_reads))
        logger.info(f"{len(tr_reads)} TR reads -> {len(centromere_reads)} centromere reads")
        return centromere_reads
    
    def hist_unit_lengths(self, tr_reads, x_min, x_max, log_scale=False):
        """Show histogram of the unit length. Counts are stratified by the total length of
        `min_ulen`-`max_ulen` bp units in a read. Only [`x_min`..`x_max`] bp units are plotted.
        """
        # Stacked histogram of unit length for each list of stratified reads
        stratified_reads = stratify_by_covered_length(tr_reads, self.min_ulen, self.max_ulen)
        traces = [make_hist([ulen for read in sreads for ulen in read_to_ulens(read, x_min, x_max)],
                            bin_size=1,
                            name=f"{covered_kb}-{covered_kb + 1}kb-covered reads")
                  for covered_kb, sreads in sorted(stratified_reads.items())]
    
        # Show peak intervals if given
        shapes = []
        if self.peak_intvls is not None:
            shapes += [make_rect(start, 0, end, 1,
                                 fill_col="gray", opacity=0.3, yref="paper", layer="above")
                       for peak_intvl in self.peak_intvls.components
                       for start, end in peak_intvl]

        layout = make_layout(title=(f"Stratified by reads according to the total length of "
                                    f"units of {self.min_ulen}-{self.max_ulen} bp"),
                             x_title="Unit length [bp]", y_title="Unit count",
                             shapes=shapes)
        layout["barmode"] = "stack"
        if log_scale:
            layout["yaxis_type"] = "log"
        show_plot(traces, layout)

    def find_peak_ulens(self, tr_reads, show_density=True):
        # Aggregate all unit lengths of [`min_ulen`..`max_ulen`] bp from reads covered
        # more than `min_covered_length` bp by such units
        all_ulens = []
        for read in tr_reads:
            ulens = read_to_ulens(read, self.min_ulen, self.max_ulen)
            if sum(ulens) >= self.min_covered_length:
                all_ulens += ulens

        # Smoothe the unit length distribution by kernel density estimation
        ulen_dens = self.smooth_distribution(all_ulens)
        if show_density:
            self.plot_density(ulen_dens)

        # Find peak unit length(s)
        self.peak_ulens = self.find_peaks(ulen_dens)
        logger.info("Peak unit lengths: " + ', '.join([f"{peak_ulen} bp"
                                                       for peak_ulen in self.peak_ulens]))

        # Compute peak unit length interval(s) from peak unit length(s)
        self.peak_intvls = interval(*[[-(- peak_ulen * (1. - self.deviation) // 1),
                                  int(peak_ulen * (1. + self.deviation))]
                                 for peak_ulen in self.peak_ulens])
        logger.info("Peak intervals: " + ', '.join([f"{start}-{end} bp"
                                                    for peak_intvl in self.peak_intvls.components
                                                    for start, end in peak_intvl]))

    def smooth_distribution(self, ulens):
        """Smooth the unit length distribution. Filtering of `unit_lens` must be finished in advance."""
        return np.exp(KernelDensity(kernel="gaussian", bandwidth=self.band_width)
                      .fit(np.array(ulens).reshape(-1, 1))
                      .score_samples(np.arange(self.min_ulen, self.max_ulen + 1).reshape(-1, 1)))

    def plot_density(self, dens):
        show_plot([make_scatter(np.arange(self.min_ulen, self.max_ulen + 1), dens,
                                mode="lines", show_legend=False)],
                  make_layout(x_title="Unit length [bp]", y_title="Density by KDE",
                              shapes=[make_line(self.min_ulen, self.min_density,
                                                self.max_ulen + 1, self.min_density,
                                                width=2, layer="above")]))

    def find_peaks(self, dens):
        """Detect peaks from the density data `dens`, a list of densities."""
        return [self.min_ulen + i for i in range(1, self.max_ulen - self.min_ulen)
                if dens[i] > dens[i - 1] and dens[i] > dens[i + 1] and dens[i] >= self.min_density]
