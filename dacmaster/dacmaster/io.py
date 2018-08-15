import os
import sys
import pickle
from logzero import logger


def load_peaks(pkl_fname):
    """
    Load a list of Peak instances.
    """

    if not os.path.isfile(pkl_fname):
        logger.error(f"{pkl_fname} does not exist. Abort.")
        sys.exit(1)

    with open(pkl_fname, 'rb') as f:
        peaks = pickle.load(f)

    if len(peaks) == 0:
        logger.error("No peak pikles were detected.")
    else:
        logger.info(f"{len(peaks)} peaks were loaded.")

    return peaks
