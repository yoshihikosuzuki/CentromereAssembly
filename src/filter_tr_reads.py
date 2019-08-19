from dataclasses import dataclass
from logzero import logger
from interval import interval
import numpy as np
from sklearn.neighbors import KernelDensity
import plotly.offline as py
import plotly.graph_objs as go
from BITS.plot.plotly import make_hist, make_scatter, make_layout, show_plot
from BITS.util.io import load_pickle
from BITS.util.interval import intvl_len


@dataclass(eq=False)
class TRReadFilter:
    """Class for filtering List[TRRead] into centromeric TR reads by finding peaks of unit lengths,
    which would be components of centromere.
    """
    tr_reads_fname : str             # output of datruf
    min_ulen       : int   = 50      # units of length in [<min_ulen>..<max_ulen>] will be used
    max_ulen       : int   = 500
    min_cover_rate : float = 0.8     # of read by all TRs
    band_width     : int   = 5       # for KDE
    min_density    : float = 0.005   # of peaks in KDE
    deviation      : float = 0.08    # <peak_ulen> * (1 +- <deviation>) will be each peak interval
    show_plot      : bool  = False

    def __post_init(self):
        self.tr_reads = load_pickle(self.tr_reads_fname)

    def run(self):
        # Find peak unit lengths from units of moderate length found in long TRs
        self.peak_ulens = self.find_peaks()
        ulens_list = ', '.join([f"{peak_ulen} bp" for peak_ulen in self.peak_ulens])
        logger.info(f"Peak unit lengths: {ulens_list}")

        # Extract reads covered by units around the peak lengths, which would come from centromere
        self.centromere_reads = self.extract_centromeric_reads()
        save_pickle(self.centromere_reads, "centromere_reads.pkl")

    def find_peaks(self):
        """Filter TR reads and units by peak unit lengths,
        based on the assumption that centromere is the major source of TRs.
        """
        # Extract reads covered by units whose lengths are within the range
        intvl = interval([self.min_ulen, self.max_ulen])
        filtered_reads = filter_reads(self.tr_reads, intvl, self.min_cover_rate)

        # Flatten lengths of all the filtered units and smooth the distribution using KDE
        ulens = [tr_unit.length for tr_read in filtered_reads for tr_unit in tr_read.units]
        ulen_dens = smooth_distribution(ulens, self.min_ulen, self.max_ulen, self.band_width)

        if self.show_plot:
            # Before smoothing
            show_plot([make_hist(ulens, start=self.min_ulen, end=self.max_ulen, bin_size=1)],
                      make_layout(None, None, x_title="Unit length", y_title="Frequency"))
            # After smoothing
            show_plot([make_scatter(np.arange(self.min_ulen, self.max_ulen + 1), ulen_dens)],
                      make_layout(None, None, x_title="Unit length", y_title="Density"))

        # Find peak unit lengths by simple sweep
        return find_peaks(ulen_dens, self.min_ulen, self.max_ulen, self.min_density)

    def extract_centromeric_reads(self):
        """Again filter the units with a range consisting of the union of each peak length +- deviation
        and filter the reads with cover rate by the filtered units.
        """
        # Redefine the acceptance range based on the peak unit lengths and 
        peak_intvls = interval(*[[-(- peak_ulen * (1. - self.deviation) // 1),
                                  int(peak_ulen * (1. + self.deviation))]
                                 for peak_ulen in self.peak_ulens])
        intvls_list = ', '.join([f"{start}-{end} bp"
                                 for peak_intvl in peak_intvls.components
                                 for start, end in peak_intvl])
        logger.info(f"Peak intervals: {intvls_list}")

        # Filter out reads with low cover rate by the units belonging to the peaks
        return filter_reads(self.tr_reads, peak_intvls, self.min_cover_rate)


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

    filtered_units = list(filter(lambda unit: unit.length in intvl, tr_read.units))
    return covered_len(filtered_units) >= min_cover_rate * tr_read.length


def smooth_distribution(X, x_min, x_max, band_width):
    """Smooth the histogram distribution made from <X> using KDE."""
    return np.exp(KernelDensity(kernel="gaussian", bandwidth=band_width)
                  .fit(X.reshape(-1, 1))
                  .score_samples(np.arange(x_min, x_max + 1).reshape(-1, 1)))


def find_peaks(dens, x_min, x_max, min_density)
    """Detect peaks from the density data <dens> = List[density] for each x in [<x_min>..<x_max>]."""
    return [x_min + i for i in range(1, x_max - x_min)
            if dens[i] > dens[i - 1] and dens[i] > dens[i + 1] and dens[i] >= min_density]
