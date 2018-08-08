import argparse
import os
import copy
import pickle
from logzero import logger

from .damaster_core import Runner


def main():
    args = load_args()

    # Check the root directory for peaks
    #if not os.path.isdir(args.peaks_dir):
    #    logger.info(f"Creating the directory: {args.peaks_dir}")
    #    run_command(f"mkdir {args.peaks_dir}")

    # Run the method of peak detection
    runner = Runner(args.unit_fasta,
                    args.peaks_dir)   # TODO: remove it
    runner.run()

    # Delte Peaks instance because no longer needed except Peak instances
    peaks = copy.deepcopy(runner.peaks)
    del runner

    # Define a set of representative monomers from somewhat homogeneous TRs
    for peak in peaks:
        peak.find_representatives()

    # Output the peaks as pickle
    # Now output the list itself instead of each peak for simplicity
    #pkl_fname = os.path.join(args.peaks_dir, "peaks.pkl")
    with open("peaks.pkl", 'wb') as f:
        pickle.dump(peaks, f)


def load_args():
    parser = argparse.ArgumentParser(
        description=("Construct master units and representative units."))

    parser.add_argument(
        "unit_fasta",
        help=("Input fasta file of the unit sequences reported by datruf. "
              "[datruf_units.fasta]"))

    parser.add_argument(
        "-n",
        "--n_core",
        type=int,
        default=1,
        help=("Degree of parallelization. [1]"))

    parser.add_argument(
        "-D",
        "--debug_mode",
        action="store_true",
        default=False,
        help=("Run in debug mode. [False]"))

    """
    parser.add_argument(   # TODO: no need of peaks dir anymore
        "-d",
        "--peaks_dir",
        default="peaks",
        help=("Directory to which dacmaster outputs results. [peaks]"))
    """

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()
