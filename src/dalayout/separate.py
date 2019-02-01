from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
import numpy as np
import pandas as pd
from BITS.run import run_edlib, run_consed, consed_to_varmat
from BITS.seq import revcomp
from BITS.utils import run_command, print_log, NoDaemonPool
import consed
from dacmaster.clustering import ClusteringVarMat


def detect_variants(encodings, repr_units, reads, peaks):
    # Collect all raw unit sequences for each (peak_id, repr_id)
    for index, df in encodings[encodings["type"] == "complete"].groupby(["peak_id", "repr_id"]):
        peak_id, repr_id = index
        fname_prefix = f"peak_{peak_id}_repr_{repr_id}.raw_units"
        with open(f"{fname_prefix}.fasta", 'w') as f:
            f.write(f">repr   length={len(repr_units[repr_id])}\n")
            f.write(f"{repr_units.loc[repr_id, 'sequennce']}\n")
            for i, e in df.iterrows():
                f.write(f">{i}/{e['strand']}/{e['start']}_{e['end']}   length={e['length']}   diff={e['diff']}\n")
                s = reads[e['read_id']]["sequence"][e["start"]:e["end"]]
                f.write(f"{s if e['strand'] == 0 else revcomp(s)}\n")

        run_command(f"grep -v '>' {fname_prefix}.fasta > {fname_prefix}.seqs")
        run_consed(f"{fname_prefix}.seqs",
                   out_prefix=f"{fname_prefix}.consed",
                   variant_vector=True,
                   variant_graph=True,
                   variant_fraction=0.0)   # TODO: parameterize
        consed_to_varmat(f"{fname_prefix}.consed")


def tsne_varmat(repr_units):
    cls = []
    for index, df in repr_units.iterrows():
        cl = ClusteringVarMat(f"peak_{df['peak_id']}_repr_{df['repr_id']}.raw_units.consed.V")
        cl.data = cl.data[:10000]   # subsampling because t-SNE is slow
        cl.N = min([cl.N, 10000])
        cl.calc_dist_mat()
        cl.plot_tsne(out_fname="peak_{df['peak_id']}_repr_{df['repr_id']}.raw_units.consed.V.tsne")
        cls.append(cl)
    return cls
