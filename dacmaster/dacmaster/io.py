import os
import sys
import pickle
from logzero import logger
from .peak import Peak
from .clustering import ClusteringSeqs


def save_peaks(peaks, pkl_fname="peaks.pkl"):
    """
    Save a list of Peak instances into a pickle file.
    """

    with open(pkl_fname, 'wb') as f:
        pickle.dump(peaks, f)


def load_peaks(pkl_fname="peaks.pkl"):
    """
    Load a list of Peak instances.
    """

    if not os.path.isfile(pkl_fname):
        logger.error(f"{pkl_fname} does not exist")
        sys.exit(1)

    with open(pkl_fname, 'rb') as f:
        peaks_loaded = pickle.load(f)

    peaks = []
    for p in peaks_loaded:
        # Peak instance
        info = p.info if hasattr(p, "info") else None
        raw_units = p.raw_units if hasattr(p, "raw_unites") else None
        peaks.append(Peak(info, raw_units))

        # instance variables of Peak class
        for attr in ("cons_units", "master_units", "master_original"):
            if hasattr(p, attr):
                setattr(peaks[-1], attr, p.attr)

        # Clustering instance
        if hasattr(peaks[-1], "cons_units"):
            setattr(peaks[-1], "cl_master", ClusteringSeqs

    if len(peaks) == 0:
        logger.warn("No peak was loaded")
    else:
        logger.info(f"{len(peaks)} peaks were loaded")
