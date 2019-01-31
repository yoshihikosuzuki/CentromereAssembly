from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
import numpy as np
import pandas as pd
from BITS.run import run_edlib
from BITS.utils import print_log, NoDaemonPool
import consed


# TODO: move
# 1. aggregation of raw units for each class from encoding DF,
# 2. calling of Consed for variant detection,
# 3. parsing of Consed output into a column of encodinf DF as variant vectors (type: global)
# (4. visualization of the units as variant vectors)

units = [set() for i in range(len(master_units))]

for read_id, encoding in encodings.items():
    read_len = len(reads[read_id])
    for start, end, length, master_id, diff, strand in encoding:
        if (start < 50) or (read_len - end < 50):   # around boundaries; maybe partial
            continue
        if diff >= 0.23:
            continue
        seq = reads[read_id][start:end]
        if strand == 1:
            seq = revcomp(seq)   # NOTE: be careful that start, end indices are those of original, forward read sequence!!!
        units[master_id].add((read_id, start, end, length, round(diff, 3), strand, seq))   # header and seq for fasta


for master_id, data in enumerate(units):
    with open(f"master_{master_id}.units.fasta", 'w') as f:
        f.write(f">master   length={len(master_units[master_id])}\n{master_units[master_id]}\n")
        for read_id, start, end, length, diff, strand, seq in sorted(data):
            f.write(f">{read_id}/{strand}/{start}_{end}   length={length}   diff={diff}\n{seq}\n")


from BITS.run import run_consed, consed_to_varmat


for master_id in range(len(master_units)):
    run_consed(f"raw_units/master_{master_id}.units.seqs", out_prefix=f"raw_units/master_{master_id}.units", variant_vector=True, variant_graph=True, variant_fraction=0.0)


for master_id in range(len(master_units)):
    consed_to_varmat(f"raw_units/master_{master_id}.units.consed")


# diff param
run_consed(f"raw_units/master_1.units.seqs", out_prefix=f"raw_units/master_1.units.t0.15", variant_vector=True, variant_graph=True, variant_fraction=0.15)
consed_to_varmat(f"raw_units/master_1.units.t0.15.consed")


from dacmaster.clustering import ClusteringVarMat


cls = []
for master_id in range(len(master_units)):
    cl = ClusteringVarMat(f"raw_units/master_{master_id}.units.consed.V")
    cl.data = cl.data[:10000]
    cl.N = min([cl.N, 10000])
    cl.calc_dist_mat()
    cl.plot_tsne()
    cls.append(cl)


# pickle save
with open("cls_varmat.pkl", 'wb') as f:
    pickle.dump(cls, f)


# num of variant sites detected for each class
!(for DATA in raw_units/*.V; do wc -l ${DATA}; done)

# After this, tried several different parameters in the jupyter notebook...
