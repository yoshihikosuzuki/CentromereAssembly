import os
from logzero import logger
import numpy as np
import pandas as pd
from scipy.spatial.distance import hamming
import plotly.graph_objs as go
from BITS.seq import homopolymer_compression
from BITS.run import run_edlib, run_consed, consed_to_varmat
from BITS.plot import generate_layout, show_plot
from BITS.utils import save_pickle, run_command
from .clustering import ClusteringSeqs, ClusteringVarMat


def _encode(read_seq, cons_seqs):
    ret = {}
    while True:
        mapping = cons_seqs.apply(lambda s: run_edlib(s,
                                                      read_seq,
                                                      "glocal",
                                                      rc=False,
                                                      rc_pruning_diff_th=0.0))
        mapping_diff = mapping.apply(lambda s: s.diff)
        if mapping_diff.min() >= 0.4:
            break
        best_idx = mapping_diff.idxmin()
        best_mapping = mapping.loc[best_idx]
        #logger.debug(f"best cons id: {best_idx}")
        #best_mapping.show()
        if best_mapping.start >= 50 and best_mapping.end <= len(read_seq) - 50:   # exclude boundary
            ret[(best_mapping.start, best_mapping.end)] = (best_idx, read_seq[best_mapping.start:best_mapping.end])
        read_seq = f"{read_seq[:best_mapping.start]}{'N' * (best_mapping.end - best_mapping.start)}{read_seq[best_mapping.end:]}"
    return [(*k, *v) for k, v in sorted(ret.items())]   # start, end, cons_id, seq


def load_graph(consed_file, panel_id, offset):
    with open(f"{consed_file}.{panel_id}.dot", 'r') as f:
        nodes, edges, variants = {}, {}, {}
        consensus = ""
        for line in f:
            line = line.strip()
            if not '[' in line:
                continue
            if not "->" in line:   # node
                node, attributes = line[:-1].split('[')
                node = node.strip()
                attributes = eval(f"dict({attributes})")
                if "label" in attributes:
                    nodes[node] = {"pos": offset + int(attributes["label"]) - 10, "color": attributes["color"]}
                else:
                    nodes[node] = {"color": attributes["color"]}
            else:   # edge
                edge, attributes = line[:-1].split('[')
                source, target = [x.strip() for x in edge.split("->")]
                attributes = eval(f"dict({attributes})")
                base, count = attributes["label"].strip().split(':')
                count = int(count)
                if nodes[source]["color"] == "blue2":   # outside the current panel
                    continue
                elif attributes["color"] == "black" and nodes[source]["color"] == "black":   # consensus path
                    if base != '-':
                        consensus += base
                elif '(' in base:   # variant path
                    base, var_id = base[:-1].split('(')
                    start = nodes[target]["pos"] if nodes[target]["color"] == "black" else None
                    end = nodes[source]["pos"] if nodes[source]["color"] == "black" else None
                    variants[(panel_id, var_id)] = (start, end, base)
    for (panel_id, var_id), (start, end, base) in variants.items():
        if start is None and end is None:
            cons = None
        elif start is None:
            cons = f"{consensus[:end]}]{consensus[end:]}"
        elif end is None:
            cons = f"{consensus[:start]}[{consensus[start:]}"
        else:
            cons = f"{consensus[:start]}[{consensus[start:end]}]{consensus[end:]}"
        variants[(panel_id, var_id)] = (start, end, base)#, cons)
    return consensus, variants


def parse_consed(consed_file, varvec=False):
    ret = {}
    ret_part = []
    with open(f"{consed_file}", 'r') as f:
        while True:
            consensus = f.readline().strip()
            if ' ' not in consensus:
                break
        ret["consensus"] = consensus
        if varvec:
            variants = []
            varmat = []
            for line in f:
                if ',' not in line:
                    continue
                a, b = line.strip().split(':')
                panel_id, var_id = map(lambda x: x.strip(), a.split(','))
                fraction, varvec = b.split()
                fraction = float(fraction)
                variants.append((panel_id, var_id, fraction))
                varmat.append(list(map(int, varvec)))
            varmat = np.array(varmat, dtype=int).T
            logger.debug(f"first varvec: {varmat[0]}")
            logger.debug(f"last varvec: {varmat[-1]}")
            varmat = varmat[1:-1]   # TODO: always exclude seed and last sequence?
            ret["varmat"] = varmat
            offset = 0
            variants_all = {}
            for panel_id in sorted(set([x[0] for x in variants])):
                cons_part, variants_part = load_graph(consed_file, panel_id, offset)
                offset += len(cons_part)
                variants_all.update(variants_part)
            #ret["ret_part"] = ret_part
            #print(variants_all)
            #print(variants)
            variants = [((panel_id, var_id), *variants_all[(panel_id, var_id)], fraction) for panel_id, var_id, fraction in variants]
            ret["variants"] = variants
    return ret


def encode(read_id, units, read_seq, hc=False, plot=False):
    logger.debug(f"{len(units)} original units")
    c = ClusteringSeqs(pd.Series([read_seq[start:end] for start, end in units]),
                       names=[f"{start}-{end}" for start, end in units],
                       cyclic=True,
                       rc=False)
    c.calc_dist_mat()
    if plot:
        c.plot_dist_mat(title="global sequence similarity")
        #c.plot_tsne()
    # Remove units which came from minor, different source
    c.cluster_hierarchical()
    c.generate_consensus()
    if c.cons_seqs.shape[0] == 0:
        logger.debug(f"No consensus units")
        return (None, None, None)
    encodings = _encode(read_seq, c.cons_seqs["sequence"])
    units = []   # TODO: XXX: name collision!!!
    variants = []
    for i, cons_seq in enumerate(c.cons_seqs["sequence"]):
        encodings_sub = [x for x in encodings if x[2] == i]
        seqs = [cons_seq] + [x[3] for x in encodings_sub]
        prefix = f"{read_id}.{i}"
        with open(f"{prefix}.seqs", 'w') as f:
            f.write('\n'.join([homopolymer_compression(s) for s in seqs] if hc else seqs) + '\n')
        run_consed(f"{prefix}.seqs", f"{prefix}.consed", variant_vector=True, variant_graph=True, variant_fraction=0.3)
        consed_to_varmat(f"{prefix}.consed")
        logger.debug(f"cons unit #{i}: {len(encodings_sub)} encoded units")
        if os.path.getsize(f"{prefix}.consed.V") == 0:   # no variants
            logger.debug(f"cons unit #{i}: no variants")
            for index, (start, end, cons_id, seq) in enumerate(encodings_sub):
                units.append((start, end, cons_id, np.empty(0)))
            variants.append([])
        else:
            logger.debug(run_command(f"cat {prefix}.consed"))
            out = parse_consed(f"{prefix}.consed", varvec=True)
            consensus = out["consensus"]
            variants.append(out["variants"])
            for (panel_id, var_id), start, end, base, fraction in out["variants"]:
                if start is None and end is None:
                    cons = None
                elif start is None:
                    cons = f"{consensus[:end]}]->[{base.upper()}]{consensus[end:]}"
                elif end is None:
                    cons = f"{consensus[:start]}[{base.upper()}]<-[{consensus[start:]}"
                else:
                    cons = f"{consensus[:start]}[{consensus[start:end]}]->[{base.upper()}]{consensus[end:]}"
                print(f"{panel_id},{var_id}")
                print(cons)
            cc = ClusteringVarMat(f"{prefix}.consed.V", [f"{x[0]}-{x[1]}" for x in encodings_sub])
            #logger.debug(f"cons unit #{i}: {cc.data.shape[1]} variants")
            for index, (start, end, cons_id, seq) in enumerate(encodings_sub):
                units.append((start, end, cons_id, cc.data[index]))
            if plot:
                cc.calc_dist_mat()
                cc.plot_dist_mat(title="units correlation")
                #cc.plot_tsne(title="units")
                # Compute correlation among variants
                cc.data = cc.data.T
                cc.calc_dist_mat()
                cc.plot_dist_mat(title="variants correlation")
                #cc.plot_tsne(title="variants")
        #run_command(f"rm {prefix}.*")
    units = sorted(units)
    logger.debug([x[:3] for x in units])
    if plot:
        s_dist_mat = [[hamming(data_i[-1], data_j[-1]) if data_i[2] == data_j[2] else 1
                       for j, data_j in enumerate(units)]
                      for i, data_i in enumerate(units)]
        trace = go.Heatmap(z=s_dist_mat,
                           colorscale="YlGnBu",
                           zmin=0,
                           zmax=1,
                           showscale=False)
        layout = generate_layout(500, 500, title="variant vector similarity")
        layout["yaxis"] = dict(autorange="reversed")
        show_plot([trace], layout)
    return (list(c.cons_seqs["sequence"]), variants, units)


def encode_read(read_id, trs, tr_reads, hc=False, plot=False):
    logger.debug(f"read {read_id}")
    return encode(read_id,
                  [unit
                   for units in trs[trs["read_id"] == read_id]["units"]
                   for unit in eval(units)],
                  tr_reads.loc[read_id]["sequence"],
                  hc=hc,
                  plot=plot)


def encode_reads(tr_file,   # output of datruf
                 tr_reads_file,   # filtered by peak unit length
                 min_n_units=50,   # in a read
                 hc=False,
                 out_fname="encoded_reads.pkl",
                 plot=False):
    trs = pd.read_csv(tr_file, sep='\t', index_col=0)
    tr_reads = pd.read_csv(tr_reads_file, sep='\t', index_col=0)

    encoded_reads = {read_id: encode_read(read_id, trs, tr_reads, hc=hc, plot=plot)
                     for peak_id, df in tr_reads.groupby("peak_id")
                     for read_id in df.index.values
                     if sum(trs[trs["read_id"] == read_id]["units"].apply(lambda s: len(eval(s)))) >= min_n_units}
    save_pickle(encoded_reads, out_fname)
