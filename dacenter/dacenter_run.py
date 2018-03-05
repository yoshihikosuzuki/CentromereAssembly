import argparse
import os
import re
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.neighbors import KernelDensity
from sklearn import cluster
from sklearn import mixture
from sklearn import decomposition
from interval import interval
import edlib

from datruf_utils import run_command

import matplotlib.pyplot as plt
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go

plt.style.use('ggplot')
#py.init_notebook_mode()   # TODO: this makes pickle in multiprocessing failed. try not to use or separate plotter


def load_units(fasta_fname):
    ret = {}
    header = ""
    seq = ""
    with open(fasta_fname, 'r') as f:
        flag_first = True
        for line in f:
            line = line.strip()
            if line[0] == '>':
                if flag_first:
                    flag_first = False
                else:
                    ret[header] = seq.upper()

                header = line[1:]
                if header in ret:
                    print("[CAUTION] Header duplication!")
                seq = ""
            else:
                seq += line

        ret[header] = seq.upper()

    return list(ret.items())


def estimate_peaks(units, min_len=50, max_len=450, band_width=10, visualize=False):
    # TODO: how to determine min/max_len? peak + density threshold?

    unit_lens = np.array([len(x[1]) for x in units])

    if visualize:
        data = [go.Histogram(x=unit_lens,
                             xbins=dict(start=20, size=1))]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Frequency"))
        py.iplot(go.Figure(data=data, layout=layout))

    ul = unit_lens[(min_len < unit_lens) & (unit_lens < max_len)]
    ls = np.linspace(min_len, max_len, max_len - min_len + 1, dtype=int)

    kde = KernelDensity(kernel='gaussian', bandwidth=band_width).fit(ul.reshape(-1, 1))
    dens = np.exp(kde.score_samples(ls.reshape(-1, 1)))

    if visualize:
        data = [go.Scatter(x=ls, y=dens)]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Density"))
        py.iplot(go.Figure(data=data, layout=layout))

    # Naive peak detection
    peaks = []
    for i in range(1, len(dens) - 1):
        if (dens[i - 1] < dens[i]) and (dens[i] > dens[i + 1]):
            peaks.append(ls[i])
            print("[INFO] peak detected:", ls[i], dens[i])

    return peaks


def calculate_peak_intervals(peaks, deviation=0.03):
    # TODO: change deviation to percent of confidence interval

    # TODO: merge multiple neighboring peaks

    peak_intervals = interval()   # from these intervals on estimated unit length, collect all unit sequences
    for x in peaks:
        peak_intervals |= interval[-(- x * (1. - deviation) // 1),
                                   int(x * (1. + deviation))]

    return peak_intervals


def output_peak_units(units, peak_intervals, out_dir):
    run_command("mkdir -p %s" % out_dir)

    writers = []
    writers_dup = []
    for i in range(1, len(peak_intervals) + 1):
        writers.append(open("%s/peak.%d.fasta" % (out_dir, i), 'w'))
        writers_dup.append(open("%s/peak.dup.%d.fasta" % (out_dir, i), 'w'))

    for header, seq in units:
        for i, intvl in enumerate(peak_intervals.components):
            if len(seq) in intvl:
                writers[i].write(">%s\n%s\n" % (header, seq))
                writers_dup[i].write(">%s\n%s\n" % (header.split('_')[0] + "_" + str(len(seq + seq)), seq + seq))
                break

    for i in range(len(writers)):
        writers[i].close()
        writers_dup[i].close()


def map_unit_to_dups(unit_seq, dup_seqs, min_dist=0.28):
    RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

    def align_edlib(a, b):
        ret = edlib.align(a, b, mode="HW", task="path")
        alignment_len = sum(list(map(int, re.findall(r'([0-9]+)', ret["cigar"]))))
        dist = ret["editDistance"] / alignment_len
        start, end = ret["locations"][0]
        end += 1
        return (dist, start, end)

    dist_array = np.full(len(dup_seqs), 0.75)
    out_seqs = []

    # TODO: store strand information for assembly

    for i, dup_seq in enumerate(dup_seqs):
        f_dist, f_start, f_end = align_edlib(unit_seq, dup_seq)

        # TODO: heuristics of not calculating revcomp if forward_dist is <20% or so

        rc_seq = "".join([RC_MAP[c] for c in dup_seq[::-1]])
        rc_dist, rc_start, rc_end = align_edlib(unit_seq, rc_seq)

        dist = min(f_dist, rc_dist)
        dist_array[i] = dist

        if dist > min_dist:   # Consed failed if units are not aligned to seed
            continue

        if f_dist <= rc_dist:
            out_seqs.append(dup_seq[f_start:f_end])
        else:
            out_seqs.append(rc_seq[rc_start:rc_end])

    return (dist_array, out_seqs)


def align_start_positions(peak_idx, peaks_dir):
    peak_units = load_units("%s/peak.%d.fasta" % (peaks_dir, peak_idx))
    dup_units = load_units("%s/peak.dup.%d.fasta" % (peaks_dir, peak_idx))
    dup_seqs = [x[1] for x in dup_units]

    # TODO: select (consensus) unit(s) that have proper boundaries

    for j, peak_unit in enumerate(peak_units):
        if j != 1:   # filter for unit number
            continue

        out_fname = "%s/peak.%d.aligned.%d.seqs" % (peaks_dir, peak_idx, j)
        if os.path.isfile(out_fname):
            print("[WARNING] %s already exists. Skip." % out_fname)
            continue

        header, seq = peak_unit
        dist_array, aligned_units = map_unit_to_dups(seq, dup_seqs)

        with open(out_fname, 'w') as f:
            for unit in aligned_units:
                f.write(unit + "\n")

    # TODO: above codes are temporary. We indeed want to output all units with
    # proper boundaries


def run_consed(units_fname, out_fname):
    command = "consed -V -w1000000 %s | awk 'NR > 5 {print $NF}'" % units_fname
    with open(out_fname, 'w')as f:
        f.write(run_command(command))


def load_vmatrix(vmatrix_fname, visualize=False):
    n_variant = int(run_command("cat %s | wc -l" % vmatrix_fname))
    n_unit = int(run_command("awk 'NR == 1 {print length($1)}' %s" % vmatrix_fname))   # NOTE: last column of vmatrix is always 0 (i.e. consensus seuqence?).

    vmatrix = np.zeros((n_variant, n_unit), dtype=int)
    with open(vmatrix_fname, 'r') as f:
        for i, line in enumerate(f):
            vmatrix[i, :] = list(map(int, list(line.strip())))

    vmatrix = vmatrix.T   # TODO: which form is better?

    if visualize:
        # Split the matrix because it is too large for single plot
        n_part = -(-n_unit // 1000)
        for i in range(n_part):
            plt.subplots(figsize=(18, 15))
            sns.heatmap(vmatrix[i * 1000:(i + 1) * 1000, :], cbar=False)
            plt.show()

            """   # plotly
            data = [go.Heatmap(z=vmatrix,
                               transpose=False,
                               colorscale='Greys',
                               reversescale=True)]
            layout = go.Layout(width=950,
                               height=950,
                               yaxis=dict(autorange='reversed'),
                               hovermode='closest')
            py.iplot(go.Figure(data=data, layout=layout))
            """

    return vmatrix


def heatmap_partition(data, cluster_index, partition_width=3):
    d = np.insert(data, 0, -1, axis=1)   # dummy data for partition lines
    order = [0] * partition_width
    for idx in set(cluster_index):
        order += np.where(cluster_index == idx)[0].tolist()
        order += [0] * partition_width
    plt.subplots(figsize=(18, 15))
    sns.heatmap(d[:, order], cbar=False)
    plt.show()


def cluster_units(vmatrix, visualize=False):
    # Hierarchical, Ward, 0.7
    hc_ward = linkage(pdist(vmatrix, metric='hamming'), method='ward')
    cluster_index = fcluster(hc_ward, 0.7, criterion='distance')

    # Another algorithm below


    if visualize:
        # Do not draw dendrogram with huge dataset
        #plt.figure(figsize=(18, 10))
        #dendrogram(hc_ward)
        #plt.show()

        if len(vmatrix[:, 0]) < 1000:   # small matrix
            heatmap_partition(vmatrix.T, cluster_index)
        else:   # huge matrix
            plt.hist(cluster_index, bins=len(set(cluster_index)))
            plt.show()
            for i in range(1, max(cluster_index) + 1):
                sns.heatmap(vmatrix[np.where(cluster_index == i)[0], :], cbar=False)
                plt.show()

    return cluster_index


def calculate_representatives(vmatrix, cluster_index, visualize=False):
    consensus = []

    for i in range(1, max(cluster_index) + 1):
        units = vmatrix[np.where(cluster_index == i)[0], :]
        cluster_size = len(units[:, 0])
        units_sum = np.sum(units, axis=0)
        consensus.append(np.array([1 if x >= cluster_size / 2 else 0 for x in units_sum]))

    if visualize:
        hc_result = linkage(pdist(consensus, metric='hamming'), method='ward')
        plt.figure(figsize=(18, 10))
        dendrogram(hc_result)
        plt.show()

    return consensus


def main():
    args = load_args()

    if args.peaks_dir is None:
        args.peaks_dir = "peaks"
        if os.path.isdir(args.peaks_dir):
            print("[ERROR] peak directory already exists.")
            exit(1)
        units = load_units(args.unit_fasta_fname)
        peaks = estimate_peaks(units)
        peak_intervals = calculate_peak_intervals(peaks)
        output_peak_units(units, peak_intervals, args.peaks_dir)
        n_peaks = len(peaks)
    else:
        if not os.path.isdir(args.peaks_dir):
            print("[ERROR] peak directory does not exist.")
            exit(1)
        command = "ls -l %s/peak.[0-9*].fasta | wc -l" % args.peaks_dir
        n_peaks = int(run_command(command))

    if not args.peaks_aligned:
        for i in range(1, n_peaks + 1):
            if i != 3:   # filter for peak number
                continue
            align_start_positions(i, args.peaks_dir)   # TODO: change to output single seqs file with units

            # NOTE: if Consed failed with a message "Cannot align to first sequence", do "$ sed -i -e "<read_id+1>d" <seqs_file>".
            vmatrix_fname = "variant_matrix.%d" % i

            run_consed(aligned_units_fname, vmatrix_fname)
            vmatrix = load_vmatrix(vmatrix_fname)
            cluster_index = cluster_units(vmatrix)
            r_units = calculate_representatives(vmatrix, cluster_index)


def load_args():
    parser = argparse.ArgumentParser(description="Run dacenter. You must "
                                     "specify one of the following arguments: "
                                     "[--unit_fasta_fname] or [--peaks_dir].")

    parser.add_argument("--unit_fasta_fname",
                        type=str,
                        default=None,
                        help=("Input fasta of unit sequences"))

    parser.add_argument("--peaks_dir",
                        type=str,
                        default=None,
                        help=("Directory in which fasta of peaks exist. "
                              "If this is specified, dacenter will not "
                              "find peaks[None]"))

    parser.add_argument("--peaks_aligned",
                        action="store_true",
                        default=False,
                        help=("Fasta of peaks are already start position"
                              "-aligned. [False]"))

    #parser.add_argument("--n_core",
    #                    type=int,
    #                    default=1,
    #                    help=("Degree of parallelization. [1]"))

    args = parser.parse_args()

    if args.unit_fasta_fname is None and args.peaks_dir is None:
        print("[ERROR] You must specify one of the following arguments:\n"
              "[--unit_fasta_fname] or [--peaks_dir].")
        exit(1)

    return args


if __name__ == "__main__":
    main()
