from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
from interval import interval
import numpy as np
import pandas as pd
from BITS.run import run_edlib
from BITS.utils import print_log, NoDaemonPool
import consed


def encode_reads(reads, repr_units):   # TODO: parallelize
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
            best_unit, best_mapping = repr_units[best_idx], mapping[best_idx]
            encodings[index] = (read_id,
                                best_mapping.start,
                                best_mapping.end,
                                best_mapping.end - best_mapping.start,
                                best_unit["peak_id"],
                                best_unit["repr_id"],
                                np.nan,
                                round(best_mapping.diff, 4),
                                best_mapping.strand)
            index += 1

            # Mask the best mapped region from the read
            read_seq = read_seq[:best_mapping.start] \
                       + "N" * (best_mapping.end - best_mapping.start) \
                       + read_seq[best_mapping.end:]
            
            # Update best mapping for each representative unit whose map region is overlapping to
            # the best mapping region at this round just above
            mapping = repr_units.apply(lambda df: mapping[df.name]
                                       if (interval([mapping[df.name].start, mapping[df.name].end])
                                           & interval([best_mapping.start, best_mapping.end]) == interval())
                                       else run_edlib(df["sequence"],
                                                      read_seq,
                                                      "glocal",
                                                      rc=True,
                                                      rc_pruning_diff_th=0.0,
                                                      strand_prior=best_mapping.strand))

    encodings = pd.from_dict(encodings,
                             orient="index",
                             columns=("read_id",
                                      "start",
                                      "end",
                                      "length",
                                      "peak_id",
                                      "repr_id",
                                      "varset_id",
                                      "identity",
                                      "strand")) \
                  .sort_values(by="start") \
                  .sort_values(by="read_id", kind="mergesort")

    return encodings
