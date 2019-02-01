from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
from interval import interval
import numpy as np
import pandas as pd
from multiprocessing import Pool
from BITS.run import run_edlib, run_consed, consed_to_varmat
from BITS.utils import print_log, NoDaemonPool, run_command
from BITS.seq import revcomp
import consed
from dacmaster.clustering import ClusteringVarMat


def encode_reads_parallel(reads, repr_units, peaks, n_core):
    unit_n_read = -(-reads.shape[0] // n_core)
    with Pool(n_core) as pool:
        rets = [ret for ret in pool.starmap(encode_reads,
                                            [(reads[i * unit_n_read:(i + 1) * unit_n_read],
                                              repr_units,
                                              peaks)
                                             for i in range(n_core)])]

    return pd.concat(rets).sort_values(by="start") \
                          .sort_values(by="read_id", kind="mergesort") \
                          .reset_index(drop=True)


def encode_reads(reads, repr_units, peaks):
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
                                               or read["length"] - best_mapping.end < 50)   # TODO: calculate from unit length
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


def cut_unit_from_read(reads, encoding):
    """
    Cut the unit sequence from a read, considering the strand.
    """

    # NOTE: <encoding> must be a single line in <encodings>
    s = reads.loc[encoding["read_id"]]["sequence"][encoding["start"]:encoding["end"]]
    return s if encoding["strand"] == 0 else revcomp(s)


def _detect_variants(peak_id, repr_id, df, repr_unit, reads):
    units = df.apply(lambda d: cut_unit_from_read(reads, d), axis=1)

    # Create a file for Consed input
    out_fname = f"peak_{peak_id}_repr_{repr_id}.raw_units"
    with open(out_fname, 'w') as f:
        f.write(repr_unit + '\n' + '\n'.join(units) + '\n')

    # Run Consed
    run_consed(out_fname,
               out_prefix=out_fname,
               variant_vector=True,
               variant_graph=True,
               variant_fraction=0.0)   # TODO: parameterize
    consed_to_varmat(f"{out_fname}.consed")   # TODO: internally recieve the matrix (implement in consed_python)

    # Load variant vectors
    with open(f"{out_fname}.consed.V", 'r') as f:
        # Use index same as <encodings> so that the results can be placed at proper locations in it
        return pd.Series(list(zip(*[list(map(int,
                                             list(line.strip()[1:-1])))
                                    for line in f])),
                         index=df.index) \
                 .apply(lambda s: np.array(s))
        # TODO: XXX: CHECK IF [1:-1] IS TRULY CORRESPOING SEQUENTIALL TO THE UNITS


def detect_variants(encodings, repr_units, reads, peaks):
    # TODO: move "peak_id" and "repr_id" in <repr_units> to index of it
    encodings["var_vec"] = pd.concat(
        [_detect_variants(peak_id,
                          repr_id,
                          df,
                          repr_units.pipe(lambda df: df[df["peak_id"] == peak_id]) \
                          .pipe(lambda df: df[df["repr_id"] == repr_id])["sequence"] \
                          .iloc[0],   # TODO: change the index of <repr_units>
                          reads)
         for (peak_id, repr_id), df
         in encodings[encodings["type"] == "complete"].groupby(["peak_id", "repr_id"])])
