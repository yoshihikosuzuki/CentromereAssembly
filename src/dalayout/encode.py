from logzero import logger
from interval import interval
import numpy as np
import pandas as pd
from multiprocessing import Pool
from BITS.run import run_edlib, run_consed, consed_to_varmat
from BITS.utils import print_log, run_command
from BITS.seq import revcomp, homopolymer_compression
import consed


def _encode_reads(repr_units, reads, peaks):
    encodings = {}
    cover_rate = {}
    index = 0
    for read_id, read in reads.iterrows():
        # Copy to a new variable so that we can modify the sequence
        read_seq = read["sequence"]
        read_len = read["length"]

        mapped_len = 0

        # Mapping of each representative unit to the read
        mapping = repr_units.apply(lambda df: run_edlib(df["sequence"],
                                                        read_seq,
                                                        "glocal",
                                                        rc=True,
                                                        rc_pruning_diff_th=0.0),
                                   axis=1)

        while True:
            mapping_diff = mapping.apply(lambda s: s.diff)

            if mapping_diff.min() >= 0.4:
                # No more good mappings
                break

            # Choose a representative unit which has the best identity to some region of the read
            best_idx = mapping_diff.idxmin()
            best_peak_id, best_repr_id = best_idx
            best_mapping = mapping.loc[best_idx]
            best_mapping_len = best_mapping.end - best_mapping.start
            th_len_from_boundary = peaks[best_peak_id].info.lens[-1] * 0.15   # TODO: search for optimal
            encodings[index] = (read_id,
                                best_mapping.start,
                                best_mapping.end,
                                best_mapping_len,
                                best_peak_id,
                                best_repr_id,
                                0,   # or "global"
                                round(best_mapping.diff, 3),
                                best_mapping.strand,
                                "boundary" if (best_mapping.start < th_len_from_boundary
                                               or read["length"] - best_mapping.end < th_len_from_boundary)
                                else "noisy" if best_mapping.diff >= 0.23
                                else "deviating" if (best_mapping_len < peaks[best_peak_id].info.min_len
                                                     or best_mapping_len > peaks[best_peak_id].info.max_len)
                                else "complete")
            index += 1
            mapped_len += best_mapping_len

            # Mask the best mapped region from the read
            read_seq = read_seq[:best_mapping.start] \
                       + "N" * (best_mapping_len) \
                       + read_seq[best_mapping.end:]

            # Update best mapping for each representative unit whose map region is overlapping to
            # the best mapping region at this round just above
            mapping = repr_units.apply(lambda df: mapping.loc[df.name]
                                       if (mapping.loc[df.name].diff >= 0.4
                                           or (interval([mapping.loc[df.name].start, mapping.loc[df.name].end])
                                               & interval([best_mapping.start, best_mapping.end]) == interval()))
                                       else run_edlib(df["sequence"],
                                                      read_seq,
                                                      "glocal",
                                                      rc=True,
                                                      rc_pruning_diff_th=0.0,
                                                      strand_prior=best_mapping.strand),
                                       axis=1)

        if mapped_len > 0:
            cover_rate[read_id] = (read_len, mapped_len, round(float(mapped_len) / read_len, 2))

    encodings = pd.DataFrame.from_dict(encodings,
                                       orient="index",
                                       columns=("read_id",
                                                "start",
                                                "end",
                                                "length",
                                                "peak_id",
                                                "repr_id",
                                                "varset_id",
                                                "diff",
                                                "strand",
                                                "type")) \
                            .sort_values(by="start") \
                            .sort_values(by="read_id", kind="mergesort") \
                            .reset_index(drop=True)

    cover_rate = pd.DataFrame.from_dict(cover_rate,
                                        orient="index",
                                        columns=("read_len",
                                                 "mapped_len",
                                                 "cover_rate"))

    return (encodings, cover_rate)


@print_log("read encoding", show_args=False)
def encode_reads(repr_units, reads, peaks, n_core):
    """
    Map the <repr_units> onto <reads> and greedily select the best-mapped unit, iteratively
    until no more good mappings are obtained.
    <reads> (pd.DataFrame) should be TR reads to reduce the computation time.
    <peaks> must be List[Peak], not PeaksFinder.
    """

    unit_n_read = -(-reads.shape[0] // n_core)
    with Pool(n_core) as pool:
        rets = [ret for ret in pool.starmap(_encode_reads,
                                            [(repr_units,
                                              reads[i * unit_n_read:(i + 1) * unit_n_read],
                                              peaks)
                                             for i in range(n_core)])]

    return (pd.concat([r[0] for r in rets]).sort_values(by="start") \
            .sort_values(by="read_id", kind="mergesort") \
            .reset_index(drop=True),
            pd.concat([r[1] for r in rets]))


def cut_unit_from_read(reads, encoding, hc):
    """
    Cut the unit sequence from a read, considering the strand.
    <encoding> must be a single line in <encodings>.
    """

    s = reads.loc[encoding["read_id"]]["sequence"][encoding["start"]:encoding["end"]]
    if hc:
        s = homopolymer_compression(s)
    return s if encoding["strand"] == 0 else revcomp(s)


def _detect_variants(peak_id, repr_id, repr_unit, reads, encoding_df, variant_fraction, hc):
    units_fname = f"consed_out/peak_{peak_id}_repr_{repr_id}.raw_units{'.hc' if hc else ''}"
    consed_fname = f"{units_fname}.t{variant_fraction}.consed"
    varmat_fname = f"{consed_fname}.V"

    # Write the representative unit and raw units aligned to it for a Consed input
    units = encoding_df.apply(lambda d: cut_unit_from_read(reads, d, hc), axis=1)
    with open(units_fname, 'w') as f:
        f.write(repr_unit + '\n' + '\n'.join(units) + '\n')

    # Run Consed
    run_consed(units_fname,
               out_fname=consed_fname,
               variant_vector=True,
               variant_graph=True,
               variant_fraction=variant_fraction)
    consed_to_varmat(consed_fname)
    # TODO: internally recieve the matrix (implement in consed_python)

    # Load variant vectors
    with open(varmat_fname, 'r') as f:
        # Use index same as <encodings> so that the results can be placed at proper locations in it
        return pd.Series(list(zip(*[list(map(int,
                                             list(line.strip()[1:-1])))
                                    for line in f])),
                         index=encoding_df.index) \
                 .apply(lambda s: np.array(s))
        # TODO: XXX: CHECK IF [1:-1] IS TRULY CORRESPOING SEQUENTIALLY TO THE UNITS BY LOOKING CONSED'S CODE


@print_log("variant detection", show_args=False)
def detect_variants(repr_units, reads, encodings, variant_fraction, hc):
    """
    Add a column of "unit variant vector" in <encodings>.
    The column is calculated for each grouped part of <encodings>, and in the end merged.
    """

    run_command("mkdir -p consed_out")
    var_vecs = [_detect_variants(peak_id,
                                 repr_id,
                                 repr_units.loc[(peak_id, repr_id), "sequence"],
                                 reads,
                                 encoding_df,
                                 variant_fraction,
                                 hc)
                for (peak_id, repr_id), encoding_df
                in encodings[encodings["type"] == "complete"].groupby(["peak_id", "repr_id"])]

    # NOTE: each layer on variant calling has a distinct column
    col_name = f"var_vec_global{variant_fraction}{'_hc' if hc else ''}"
    encodings[col_name] = pd.concat(var_vecs)
