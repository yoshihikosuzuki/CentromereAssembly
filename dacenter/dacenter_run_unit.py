import argparse
import os

from dacenter_unit import Clustering
from datruf_utils import run_command


def run_consed(in_fname):
    out_fname_prefix = os.path.splitext(in_fname)[0]

    print("[INFO] If Consed failed with \"Cannot align to first sequence\", "
          "try \"$ sed -i -e \"<read_id + 1>d\" <aligned_units_fname>\".")

    # Consed output with variant matrix
    command = "consed -V -w1000000 %s" % in_fname
    consed_fname = out_fname_prefix + ".consed"
    with open(consed_fname, 'w')as f:
        f.write(run_command(command))

    # Extract consensus sequence
    command = ("awk 'BEGIN {seq = \"\"} $0 == \"\" {exit} {seq = seq $0} END {print seq}' %s"
               % consed_fname)
    cons_seq = run_command(command).strip()
    with open(out_fname_prefix + ".consensus.fasta", 'w') as f:
        f.write(">consensus/0/0_%d\n%s\n" % (len(cons_seq), cons_seq))

    # Extract variant matrix
    command = "awk 'NR > 5 {print $NF}' %s" % consed_fname
    with open(out_fname_prefix + ".V", 'w') as f:
        f.write(run_command(command))

    return out_fname_prefix + ".V"


def main():
    args = load_args()

    # TODO: generate variant graph as well?
    if args.aligned_units_fname is not None:
        args.variant_matrix_fname = run_consed(args.aligned_units_fname)

    clustering = Clustering(args.variant_matrix_fname)
    clustering.ward()
    print("[INFO] %d clusters were generated." % len(set(clustering.assignment)))
    clustering.generate_representative_units()
    clustering.output_representative_units("repr.ward.fasta")


def load_args():
    parser = argparse.ArgumentParser(
        description=("Run dacenter for clustering of units. You must specify "
                     "one of the following arguments: [--aligned_units_fname] "
                     "or [--variant_matrix_fname]."))

    parser.add_argument(
        "--aligned_units_fname",
        default=None,
        help=(".seq file of the (start-aligned) unit sequences [None]"))

    parser.add_argument(
        "--variant_matrix_fname",
        default=None,
        help=("Precomputed variant matrix file [None]"))

    args = parser.parse_args()

    if args.aligned_units_fname is None and args.variant_matrix_fname is None:
        print("[ERROR] You must specify one of the following arguments:\n"
              "[--aligned_units_fname] or [--variant_matrix_fname].")
        exit(1)

    return args


if __name__ == "__main__":
    main()
