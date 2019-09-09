from dataclasses import dataclass
from logzero import logger
from interval import interval
import numpy as np
from sklearn.neighbors import KernelDensity
from BITS.plot.plotly import make_rect, make_hist, make_scatter, make_layout, show_plot
from BITS.util.io import load_pickle, save_pickle
from BITS.util.interval import intvl_len


@dataclass(eq=False)
class TRReadFilter:
    """Class for filtering List[TRRead] into centromeric TR reads by finding peaks of unit lengths,
    which should be a sign of centromeric TR units.

    The flow of the filtering is as follows:
      - Given: TR reads (= reads w/ TR(s) of any unit length & any copy number)
            => TR-contained reads (= reads contained within TR(s))
            => Centromeric reads (= reads contained within TR(s) and having units of peak lengths)

    Before executing run(), It is recomended to adjust the parameters through looking at histograms
    of unit length using hist_all_units() and hist_filtered_units() methods inside Jupyter Notebook.
    """
    tr_reads_fname : str                # output of datruf
    min_ulen       : int      = 50      # units of length in [<min_ulen>..<max_ulen>] will be used
    max_ulen       : int      = 500
    min_cover_rate : float    = 0.8     # of read by all TRs
    band_width     : int      = 5       # for KDE
    min_density    : float    = 0.005   # of peaks in KDE
    deviation      : float    = 0.1     # <peak_ulen> * (1 +- <deviation>) will be each peak interval

    def __post_init__(self):
        self.tr_reads = load_pickle(self.tr_reads_fname)

    def run(self):
        # Find peak lengths and intervals of TR untis from TR-contained reads
        peak_intvls = self.find_peaks()

        # Again filter the TR reads using the peak intervals
        centromere_reads = filter_reads(self.tr_reads, peak_intvls, self.min_cover_rate)
        logger.info(f"{len(self.tr_reads)} TR reads -> {len(centromere_reads)} centromere reads")
        save_pickle(centromere_reads, "centromere_reads.pkl")

    def find_peaks(self, plot=False):
        """Filter TR reads and units by peak unit lengths,
        based on the assumption that centromere is the major source of TRs.
        """
        # Extract reads covered by units whose lengths are within the range
        filtered_reads = filter_reads(self.tr_reads,
                                      interval([self.min_ulen, self.max_ulen]),
                                      self.min_cover_rate)

        # Flatten lengths of all the filtered units and smooth the distribution using KDE
        ulens = [tr_unit.length for tr_read in filtered_reads for tr_unit in tr_read.units]
        ulen_dens = smooth_distribution(ulens, self.min_ulen, self.max_ulen, self.band_width)

        # Find peak unit lengths by simple sweep
        peak_ulens = find_peaks(ulen_dens, self.min_ulen, self.max_ulen, self.min_density)
        logger.info("Peak unit lengths: " + ', '.join([f"{peak_ulen} bp" for peak_ulen in peak_ulens]))

        # Compute peak unit length intervals permitting some <deviation>
        peak_intvls = interval(*[[-(- peak_ulen * (1. - self.deviation) // 1),
                                  int(peak_ulen * (1. + self.deviation))]
                                 for peak_ulen in peak_ulens])
        logger.info("Peak intervals: " + ', '.join([f"{start}-{end} bp"
                                                    for peak_intvl in peak_intvls.components
                                                    for start, end in peak_intvl]))

        if plot:
            peak_rects = [make_rect(start, 0, end, 1, yref="paper", layer="above")
                          for peak_intvl in peak_intvls.components for start, end in peak_intvl]
            # Before smoothing
            show_plot([make_hist(ulens, start=self.min_ulen, end=self.max_ulen, bin_size=1)],
                      make_layout(x_title="Unit length [Filtered]", y_title="Frequency", shapes=peak_rects))
            # After smoothing
            show_plot([make_scatter(np.arange(self.min_ulen, self.max_ulen + 1), ulen_dens,
                                    mode="lines", show_legend=False)],
                      make_layout(x_title="Unit length [Filtered]", y_title="Density by KDE"))

        return peak_intvls

    def hist_all_units(self):
        """Show unit length histogram with all units."""
        all_ulens = [tr_unit.length for tr_read in self.tr_reads for tr_unit in tr_read.units]
        show_plot([make_hist(all_ulens, start=self.min_ulen, end=self.max_ulen, bin_size=1)],
                  make_layout(x_title="Unit length [All]", y_title="Frequency"))

    def hist_filtered_units(self):
        """Show unit length histogram with filtered units. Units inside the shades will be output."""
        self.find_peaks(plot=True)


def filter_reads(tr_reads, intvl, min_cover_rate):
    """Extract reads covered more than <min_cover_rate> by units whose lengths are within <intvl>."""
    return list(filter(lambda tr_read: is_eligible(tr_read, intvl, min_cover_rate),
                       tr_reads))


def is_eligible(tr_read, intvl, min_cover_rate):
    """Judge whether or not <tr_read> is covered by units of length within <intvl> at the percentage of
    <min_cover_rate>.
    """
    def covered_len(units):
        return intvl_len(interval(*[(unit.start, unit.end - 1) for unit in units]))

    # First filter by unit length
    filtered_units = list(filter(lambda unit: unit.length in intvl, tr_read.units))
    # Then filter by cover rate
    return covered_len(filtered_units) >= min_cover_rate * tr_read.length


def smooth_distribution(X, x_min, x_max, band_width):
    """Smooth the histogram distribution made from <X> using KDE."""
    return np.exp(KernelDensity(kernel="gaussian", bandwidth=band_width)
                  .fit(np.array(X).reshape(-1, 1))
                  .score_samples(np.arange(x_min, x_max + 1).reshape(-1, 1)))


def find_peaks(dens, x_min, x_max, min_density):
    """Detect peaks from the density data <dens> = List[density] for each x in [<x_min>..<x_max>]."""
    return [x_min + i for i in range(1, x_max - x_min)
            if dens[i] > dens[i - 1] and dens[i] > dens[i + 1] and dens[i] >= min_density]
