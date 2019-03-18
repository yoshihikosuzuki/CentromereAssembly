from logzero import logger
import numpy as np
import pandas as pd
from BITS.run import run_edlib, run_consed, consed_to_varmat
from .clustering import ClusteringSeqs, ClusteringVarMat


def _encode(read_seq, cons_seq):
    ret = {}
    while True:
        mapping = run_edlib(cons_seq, read_seq, "glocal", rc=False, rc_pruning_diff_th=0.0)
        if mapping.diff >= 0.4:
            break
        mapping.show()
        if mapping.start >= 50 and mapping.end <= len(read_seq) - 50:   # exclude boundary
            ret[(mapping.start, mapping.end)] = read_seq[mapping.start:mapping.end]
        read_seq = f"{read_seq[:mapping.start]}{'N' * (mapping.end - mapping.start)}{read_seq[mapping.end:]}"
    return [(*k, v) for k, v in sorted(ret.items())]   # start, end, seq


def encode(units, read_seq, plot=False):
    logger.debug(f"{len(units)} original units")
    c = ClusteringSeqs([read_seq[start:end] for start, end in units],
                       names=[f"{start}-{end}" for start, end in units],
                       cyclic=True,
                       rc=False)
    c.calc_dist_mat()
    if plot:
        c.plot_dist_mat()
        c.plot_tsne()
    c.generate_consensus()
    cons_seq = c.cons_seqs["sequence"][0]
    encodings = _encode(read_seq, cons_seq)
    seqs = [cons_seq] + [x[2] for x in encodings]
    with open("tmp.seqs", 'w') as f:
        f.write('\n'.join(seqs) + '\n')
    run_consed("tmp.seqs", "tmp.consed", variant_vector=True, variant_graph=True, variant_fraction=0.3)
    consed_to_varmat("tmp.consed")
    logger.debug(f"{len(encodings)} encoded units")
    c = ClusteringVarMat("tmp.consed.V", [f"{x[0]}-{x[1]}" for x in encodings])   # TODO: modify ClusteringVarMat so that it accepts custom names
    logger.debug(f"{c.data.shape[1]} variants")
    c.calc_dist_mat()
    if plot:
        c.plot_dist_mat()
        c.plot_tsne()


def encode_read(read_id, trs, tr_reads, plot=False):
    encode([unit
            for units in trs[trs["read_id"] == read_id]["units"]
            for unit in eval(units)],
           tr_reads.loc[read_id]["sequence"],
           plot=plot)


def encode_reads(tr_file,   # output of datruf
                 tr_reads_file):   # filtered by peak unit length
    trs = pd.read_csv(tr_file, sep='\t', index_col=0)
    tr_reads = pd.read_csv(tr_reads_file, sep='\t', index_col=0)

    for peak_id, df in tr_reads.groupby("peak_id"):
        for read_id in df.index.values:
            encode_read(read_id, trs, tr_reads)
