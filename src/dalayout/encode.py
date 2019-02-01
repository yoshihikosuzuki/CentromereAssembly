from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
from interval import interval
import numpy as np
import pandas as pd
from BITS.run import run_edlib, run_consed, consed_to_varmat
from BITS.utils import print_log, NoDaemonPool, run_command
from BITS.seq import revcomp
import consed
from dacmaster.clustering import ClusteringVarMat


def encode_reads(reads, repr_units, peaks):   # TODO: parallelize
    """
    Map the <repr_units> onto <reads> and greedily select the best-mapped unit, iteratively
    until no more good mappings are obtained.
    <reads> (pd.DataFrame) should be TR reads to reduce the computation time.
    <peaks> must be List[Peak], not PeaksFinder.
    """

    encodings = {}
    index = 0
    for read_id, read in reads.iterrows():
        # Copy to a new variable so that we can modify the sequence
        read_seq = read["sequence"]

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
            best_unit, best_mapping = repr_units.loc[best_idx], mapping.loc[best_idx]
            encodings[index] = (read_id,
                                best_mapping.start,
                                best_mapping.end,
                                best_mapping.end - best_mapping.start,
                                best_unit["peak_id"],
                                best_unit["repr_id"],
                                0,   # or "global"
                                round(best_mapping.diff, 3),
                                best_mapping.strand,
                                "boundary" if (best_mapping.start < 50
                                               and read["length"] - best_mapping.end < 50)   # TODO: calculate from unit length
                                else "noisy" if best_mapping.diff >= 0.23
                                else "deviating" if (best_mapping.end - best_mapping.start < peaks[best_unit["peak_id"]].info.min_len
                                                     or best_mapping.end - best_mapping.start > peaks[best_unit["peak_id"]].info.max_len)   # TODO: check if this category exists instead of "noisy"
                                else "complete")
            index += 1

            # Mask the best mapped region from the read
            read_seq = read_seq[:best_mapping.start] \
                       + "N" * (best_mapping.end - best_mapping.start) \
                       + read_seq[best_mapping.end:]

            # Update best mapping for each representative unit whose map region is overlapping to
            # the best mapping region at this round just above
            mapping = repr_units.apply(lambda df: mapping.loc[df.name]
                                       if (interval([mapping.loc[df.name].start, mapping.loc[df.name].end])
                                           & interval([best_mapping.start, best_mapping.end]) == interval())
                                       else run_edlib(df["sequence"],
                                                      read_seq,
                                                      "glocal",
                                                      rc=True,
                                                      rc_pruning_diff_th=0.0,
                                                      strand_prior=best_mapping.strand),
                                       axis=1)

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

    return encodings


def detect_variants(encodings, repr_units, reads, peaks, plot_tsne=False):
    """
    From <encodings>, collect all raw unit sequences and throw into Consed
    for each pair of (peak_id, repr_id).
    """

    for index, df in encodings[encodings["type"] == "complete"].groupby(["peak_id", "repr_id"]):
        peak_id, repr_id = index
        fname_prefix = f"peak_{peak_id}_repr_{repr_id}.raw_units"
        with open(f"{fname_prefix}.fasta", 'w') as f:
            f.write(f">repr   length={len(repr_units.loc[repr_id])}\n")
            f.write(f"{repr_units.loc[repr_id, 'sequence']}\n")
            for i, e in df.iterrows():
                f.write(f">{i}/{e['strand']}/{e['start']}_{e['end']}   length={e['length']}   diff={e['diff']}\n")
                s = reads.loc[e['read_id']]["sequence"][e["start"]:e["end"]]
                f.write(f"{s if e['strand'] == 0 else revcomp(s)}\n")

        run_command(f"grep -v '>' {fname_prefix}.fasta > {fname_prefix}.seqs")
        run_consed(f"{fname_prefix}.seqs",
                   out_prefix=fname_prefix,
                   variant_vector=True,
                   variant_graph=True,
                   variant_fraction=0.0)   # TODO: parameterize
        consed_to_varmat(f"{fname_prefix}.consed")   # TODO: extract positions information and so on

        if plot_tsne:
            cl = ClusteringVarMat(f"{fname_prefix}.consed.V")
            cl.data = cl.data[:10000]   # subsampling because t-SNE is slow
            cl.N = min([cl.N, 10000])
            cl.calc_dist_mat()
            cl.plot_tsne(out_fname="{fname_prefix}.consed.V.tsne.png")
