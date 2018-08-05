import argparse
import os
import copy
import pickle
from logzero import logger

from BITS.utils import run_command
from .dacmaster_core import Runner


def main():
    args = load_args()

    # Check the root directory for peaks
    if not os.path.isdir(args.peaks_dir):
        logger.info(f"Creating the directory: {args.peaks_dir}")
        run_command(f"mkdir {args.peaks_dir}")

    # Run the method of peak detection
    runner = Runner(args.unit_fasta,
                    args.peaks_dir)
    runner.run()

    # Delte Peaks instance because no longer needed except Peak instances
    peaks = copy.deepcopy(runner.peaks)
    del runner

    # Define a set of representative monomers from somewhat homogeneous TRs
    for peak in peaks:
        peak.find_representatives()

    # Output the peaks as pickle
    # Now output the list itself instead of each peak for simplicity
    pkl_fname = os.path.join(args.peaks_dir, "peaks.pkl")
    with open(pkl_fname, 'wb') as f:
        pickle.dump(peaks, f)


def load_args():
    parser = argparse.ArgumentParser(
        description=("Run dacmaster."))

    parser.add_argument(
        "unit_fasta",
        help=("Input fasta file of the unit sequences reported by datruf. "
              "[datruf_units.fasta]"))

    parser.add_argument(
        "-d",
        "--peaks_dir",
        default="peaks",
        help=("Directory to which dacmaster outputs results. [peaks]"))

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()
