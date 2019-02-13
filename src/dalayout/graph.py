from typing import List
from dataclasses import dataclass, field, InitVar
from multiprocessing import Pool
from logzero import logger
import numpy as np
import pandas as pd
import igraph as ig
from collections import Counter
import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
from BITS.utils import run_command, sge_nize, save_pickle, make_line

plt.style.use('ggplot')


def to_unit_type(s, strand):
    ret = list(s[["peak_id", "repr_id", "strand", "type"]])
    if strand == 1:
        ret[2] = strand - ret[2]
    return ret


def plot_alignment_mat(read_df_i, read_df_j, strand, score_mat, dp, path):
    # NOTE: default heatmap is [x (on y-axis) * y (on x-axis)] for a distance matrix [x * y],
    #       so transpose the matrix.
    
    # Score matrix
    trace1 = go.Heatmap(z=score_mat.T,
                       text=np.array([[f"{ri}: {to_unit_type(read_df_i.iloc[ri], 0)} vs {ci}: {to_unit_type(read_df_j.iloc[ci], strand)}<br>%sim={c:.2}"
                                       for ci, c in enumerate(r)]
                                      for ri, r in enumerate(score_mat)]).T,
                       hoverinfo="text",
                       colorscale='YlGnBu',
                       #colorscale='Greys',
                       zmin=0.6,
                       zmax=1,
                       reversescale=True,
                       showscale=False)

    trace2 = go.Scatter(x=[x[0] for x in path],
                        y=[x[1] for x in path],
                        text=[f"{x[0]}: {to_unit_type(read_df_i.iloc[x[0]], 0)} vs {x[1]}: {to_unit_type(read_df_j.iloc[x[1]], strand)}<br>%sim={score_mat[x[0]][x[1]]:.2}"
                              for x in path],
                        hoverinfo="text",
                        name="optimal path")

    layout = go.Layout(width=800,
                       height=800,
                       xaxis=dict(
                           showgrid=False,
                           zeroline=False
                       ),
                       yaxis=dict(
                           autorange="reversed",
                           scaleanchor="x",
                           showgrid=False,
                           zeroline=False
                       ),
                       showlegend=True)

    py.iplot(go.Figure(data=[trace1, trace2], layout=layout))
    
    # DP matrix
    trace3 = go.Heatmap(z=dp.T,
                       text=np.array([[f"{ri - 1}: {to_unit_type(read_df_i.iloc[ri - 1], 0)} vs {ci - 1}: {to_unit_type(read_df_j.iloc[ci - 1], strand)}<br>%sim={c:.2}"
                                       if ri * ci != 0 else "0"
                                       for ci, c in enumerate(r)]
                                      for ri, r in enumerate(dp)]).T,
                       hoverinfo="text",
                       colorscale='YlGnBu',
                       #colorscale='Greys',
                       reversescale=True,
                       showscale=False)

    trace4 = go.Scatter(x=[x[0] + 1 for x in path],
                        y=[x[1] + 1 for x in path],
                        text=[f"{x[0]}: {to_unit_type(read_df_i.iloc[x[0]], 0)} vs {x[1]}: {to_unit_type(read_df_j.iloc[x[1]], strand)}<br>%sim={dp[x[0] + 1][x[1] + 1]:.2}"
                              for x in path],
                        hoverinfo="text",
                        name="optimal path")

    layout2 = go.Layout(width=800,
                       height=800,
                       xaxis=dict(
                           showgrid=False,
                           zeroline=False
                       ),
                       yaxis=dict(
                           autorange="reversed",
                           scaleanchor="x",
                           showgrid=False,
                           zeroline=False
                       ),
                       showlegend=True)

    py.iplot(go.Figure(data=[trace3, trace4], layout=layout2))


def calc_score(si, sj, strand, varvec_colname, match_th):
    if si["peak_id"] != sj["peak_id"]:   # unit length mismatch
        return 0
    elif abs(si["strand"] - sj["strand"]) != strand:   # strand mismatch
        return 0
    elif si["type"] == "boundary" or sj["type"] == "boundary":
        return match_th
    elif si["type"] == "deviating" or sj["type"] == "deviating":
        if si["type"] == sj["type"]:   # both units are deviating; this would be a signature rather than error
            return 1   # TODO: consider more
        else:
            return match_th   # no positive nor negative score
    elif si["repr_id"] != sj["repr_id"]:
        return 0   # TODO: prepare score matrix based on the global similarity matrix
    elif si["type"] == "noisy" or sj["type"] == "noisy":
        return match_th - 0.1   # XXX: consider more   TODO: permit repr_id mismatch for noisy units? (though it is dangerous)
    else:   # same peak, same repr, same strand, complete unit
        assert si[varvec_colname].shape[0] == sj[varvec_colname].shape[0], "Variant vector length mismatch"
        return 1. - float(np.count_nonzero(si[varvec_colname] != sj[varvec_colname])) / si[varvec_colname].shape[0]


def calc_score_mat(read_df_i, read_df_j, strand, varvec_colname, match_th):
    return np.array([[calc_score(read_df_i.iloc[i], read_df_j.iloc[j], strand, varvec_colname, match_th)
                      for j in range(read_df_j.shape[0])]
                     for i in range(read_df_i.shape[0])])


def calc_alignment(score_mat, match_th, indel_penalty):
    # Fill DP matrix
    dp = np.zeros((score_mat.shape[0] + 1, score_mat.shape[1] + 1))
    for i in range(1, dp.shape[0]):
        for j in range(1, dp.shape[1]):
            dp[i][j] = max([dp[i - 1][j - 1] + score_mat[i - 1][j - 1] - match_th,
                            dp[i - 1][j] - indel_penalty,
                            dp[i][j - 1] - indel_penalty])   # TODO: reconsider the scoring system

    # Find the starting point of traceback
    if np.max(dp[-1][1:]) >= np.max(dp.T[-1][1:]):   # maximum is on the end row
        argmax = np.array([dp.shape[0] - 1, np.argmax(dp[-1][1:]) + 1])
    else:   # on the end column
        argmax = np.array([np.argmax(dp.T[-1][1:]) + 1, dp.shape[1] - 1])

    # Traceback
    path = [argmax]   # [(unit_i, unit_j), ...] from the traceback starting point
    while True:
        if argmax[0] == 1 or argmax[1] == 1:   # reached the start row or colum
            break
        diag = dp[argmax[0] - 1][argmax[1] - 1] + score_mat[argmax[0] - 1][argmax[1] - 1] - match_th
        horizontal = dp[argmax[0]][argmax[1] - 1] - indel_penalty
        vertical = dp[argmax[0] - 1][argmax[1]] - indel_penalty
        maximum = np.argmax([diag, horizontal, vertical])   # pointer to the next cell
        if maximum == 0:   # diagonal
            argmax = argmax - 1
        elif maximum == 1:   # horizonal
            argmax = argmax - [0, 1]   # copy object
        else:   # vertical
            argmax = argmax - [1, 0]
        path = [argmax] + path   # add the next cell at the FRONT of the list

    score = dp[path[-1][0]][path[-1][1]]
    path = [p - 1 for p in path]   # -1 for converting to the indices in read_df and score_mat
    return (dp, path, score)


def _svs_read_alignment(read_i,
                        read_j,
                        strand,
                        read_df_i,
                        read_df_j,
                        varvec_colname,
                        match_th=0.75,
                        indel_penalty=0.2,
                        plot=False):

    strand = 0 if strand == 'f' else 1   # for strand match judgement between two units

    # Calculate match score for every unit pair between read_i and read_j
    score_mat = calc_score_mat(read_df_i,
                               read_df_j if strand == 0 else read_df_j[::-1],
                               strand,
                               varvec_colname,
                               match_th)

    # Take alignment by solving a DP which intermediates between NW and SW
    dp, path, score = calc_alignment(score_mat, match_th, indel_penalty)
    mean_score = score / len(path)

    n_unit_i, n_unit_j = read_df_i.shape[0], read_df_j.shape[0]
    start_unit_i, start_unit_j = path[0]
    end_unit_i, end_unit_j = path[-1]
    if strand == 1:
        start_unit_j, end_unit_j = n_unit_j - 1 - end_unit_j, n_unit_j - 1 - start_unit_j
    start_bp_i, start_bp_j = read_df_i.iloc[start_unit_i]["start"], read_df_j.iloc[start_unit_j]["start"]
    end_bp_i, end_bp_j = read_df_i.iloc[end_unit_i]["end"], read_df_j.iloc[end_unit_j]["end"]

    alignment_bp_i = end_bp_i - start_bp_i
    alignment_bp_j = end_bp_j - start_bp_j
    overlap_len = sum([max(read_df_i.iloc[i]["length"], read_df_j.iloc[j]["length"]) for i, j in path])   # XXX: TODO: compute accurate value!!!

    if plot:
        plot_alignment_mat(read_df_i,
                           read_df_j if strand == 0 else read_df_j[::-1],
                           strand,
                           score_mat,
                           dp,
                           path)

    return [read_i,
            read_j,
            strand,
            start_unit_i,
            end_unit_i,
            n_unit_i,
            start_unit_j,
            end_unit_j,
            n_unit_j,
            start_bp_i,
            end_bp_i,
            # TODO: read_len_i
            start_bp_j,
            end_bp_j,
            # TODO: read_len_j
            score,
            mean_score,
            alignment_bp_i,
            alignment_bp_j,
            overlap_len,
            path]


def svs_read_alignment(read_i,
                       read_j,
                       strand,
                       read_df_i,
                       read_df_j,
                       comp_i,
                       comp_j,
                       varvec_colname,
                       th_n_shared_units,
                       th_mean_score,
                       th_ovlp_len):

    if sum((comp_i & comp_j).values()) >= th_n_shared_units:
        overlap = _svs_read_alignment(read_i, read_j, strand, read_df_i, read_df_j, varvec_colname)
        if overlap[11] >= th_mean_score and overlap[14] >= th_ovlp_len:
            return overlap
    return None


def svs_read_alignment_mult(list_pairs,
                            read_dfs,
                            comps,
                            varvec_colname,
                            th_n_shared_units,
                            th_mean_score,
                            th_ovlp_len):

    return [svs_read_alignment(read_i,
                               read_j,
                               strand,
                               read_dfs[read_i],
                               read_dfs[read_j],
                               comps[(read_i, 'f')],
                               comps[(read_j, strand)],
                               varvec_colname,
                               th_n_shared_units,
                               th_mean_score,
                               th_ovlp_len)
            for read_i, read_j, strand in list_pairs]


@dataclass(repr=False, eq=False)
class Overlap:
    encodings: pd.DataFrame
    varvec_colname: str = "var_vec_global0.0"
    th_n_shared_units: int = 10   # for preliminary filtering   # TODO: change to bp
    th_mean_score: float = 0.01
    th_ovlp_len: int = 1000   # in bp

    read_ids: List[int] = field(init=False)
    read_dfs: dict = field(init=False)
    comps: dict = field(init=False)

    def __post_init__(self):
        # TODO: at least we should remove reads only with noisy units?
        self.read_dfs = dict(tuple(self.encodings.groupby("read_id")))
        self.read_ids = sorted(self.read_dfs.keys())
        self.comps = {(read_id, strand): self.df_to_composition(self.read_dfs[read_id], strand)
                      for read_id in self.read_ids
                      for strand in ['f', 'r']}

    def df_to_composition(self, df, strand):
        return Counter(df.apply(lambda df: (df["peak_id"],
                                            df["repr_id"],
                                            df["strand"] if strand == 'f' else 1 - df["strand"]),
                                axis=1))

    def ava_read_alignment(self, n_core=1):
        """
        Non-distributed version of all-vs-all read overlap calculation.
        """

        # Pairs of reads to be aligned
        list_pairs = [(read_i, read_j, strand)
                      for i, read_i in enumerate(self.read_ids)
                      for read_j in self.read_ids[i + 1:]
                      for strand in ['f', 'r']]
        self._ava_read_alignment(list_pairs, n_core)

    def ava_read_alignment_distribute(self, n_distribute, n_core):
        """
        Distributed all-vs-all read overlap calculation.
        """

        list_pairs = [(read_i, read_j, strand)
                      for i, read_i in enumerate(self.read_ids)
                      for read_j in self.read_ids[i + 1:]
                      for strand in ['f', 'r']]
        n_pairs = len(list_pairs)
        n_sub = -(-n_pairs // n_distribute)
        list_pairs_sub = [list_pairs[i * n_sub:(i + 1) * n_sub - 1]
                          for i in range(n_distribute)]
        save_pickle(self, "overlap_obj.pkl")
        n_digit = int(np.log10(n_distribute) + 1)
        for i, lp in enumerate(list_pairs_sub):
            index = str(i + 1).zfill(n_digit)
            save_pickle(lp, f"list_pairs.{index}.pkl")
            script_fname = f"ava_pair.sge.{index}"
            with open(script_fname, 'w') as f:
                script = ' '.join([f"run_ava_sub.py",
                                   f"-n {n_core}",
                                   f"overlap_obj.pkl",
                                   f"list_pairs.{index}.pkl",
                                   f"overlaps.{index}.pkl"])
                script = sge_nize(script,
                                  job_name="run_ava_sub",
                                  n_core=n_core,
                                  wait=False)
                f.write(script)
            run_command(f"qsub {script_fname}")

    def _ava_read_alignment(self, list_pairs, n_core):
        n_pairs = len(list_pairs)
        n_sub = -(-n_pairs // n_core)
        list_pairs_sub = [(list_pairs[i * n_sub:(i + 1) * n_sub - 1],
                           self.read_dfs,
                           self.comps,
                           self.varvec_colname,
                           self.th_n_shared_units,
                           self.th_mean_score,
                           self.th_ovlp_len)
                          for i in range(n_core)]

        overlaps = {}
        index = 0
        with Pool(n_core) as pool:
            for ret in pool.starmap(svs_read_alignment_mult, list_pairs_sub):
                for r in ret:
                    if r is not None:
                        overlaps[index] = r
                        index += 1

        self.overlaps = pd.DataFrame.from_dict(overlaps,
                                               orient="index",
                                               columns=("read_i",
                                                        "read_j",
                                                        "strand",
                                                        "i_start_unit",
                                                        "i_end_unit",
                                                        "i_n_unit",
                                                        "j_start_unit",
                                                        "j_end_unit",
                                                        "j_n_unit",
                                                        "i_start_bp",
                                                        "i_end_bp",
                                                        "j_start_bp",
                                                        "j_end_bp",
                                                        "score",
                                                        "mean_score",
                                                        "i_alignment_bp",
                                                        "j_alignment_bp",
                                                        "ovelap_len",
                                                        "path")) \
                                    .sort_values(by="strand") \
                                    .sort_values(by="read_j", kind="mergesort") \
                                    .sort_values(by="read_i", kind="mergesort") \
                                    .reset_index(drop=True)


def construct_string_graph(overlaps, repr_units, th_mean_score=0.0, th_overlap_len=3000):
    def alignment_to_length(a):
        return sum([repr_units.loc[(y[0], y[1])]["length"] for x, y, z in a])

    # NOTE: Because igraph prefers static graph construction, first list the vertices and edges up.
    nodes, edges = set(), set()
    for i, overlap in overlaps.iterrows():
        f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l, score, mean_score, alignment = overlap
        if mean_score < th_mean_score or alignment_to_length(alignment) < th_overlap_len:
            continue

        if strand == "r":  # reversed alignment, swapping the begin and end coordinates
            g_b, g_e = g_e, g_b

        # build the string graph edges for each overlap
        if f_b > 3:
            if g_b < g_e:
                """
                     f.B         f.E
                  f  ----------->
                  g         ------------->
                            g.B           g.E
                """
                if f_b < 3 or g_l - g_e < 3:   # contained
                    continue
                nodes.update(["%s:B" % g_id,
                              "%s:B" % f_id,
                              "%s:E" % f_id,
                              "%s:E" % g_id])
                edges.update([("%s:B" % g_id, "%s:B" % f_id),
                              ("%s:E" % f_id, "%s:E" % g_id)])
            else:
                """
                     f.B         f.E
                  f  ----------->
                  g         <-------------
                            g.E           g.B
                """
                if f_b < 3 or g_e < 3:
                    continue
                nodes.update(["%s:E" % g_id,
                              "%s:B" % f_id,
                              "%s:E" % f_id,
                              "%s:B" % g_id])
                edges.update([("%s:E" % g_id, "%s:B" % f_id),
                              ("%s:E" % f_id, "%s:B" % g_id)])
        else:
            if g_b < g_e:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         ------------->
                            g.B           g.E
                """
                if g_b < 3 or f_l - f_e < 3:
                    continue
                nodes.update(["%s:B" % f_id,
                              "%s:B" % g_id,
                              "%s:E" % g_id,
                              "%s:E" % f_id])
                edges.update([("%s:B" % f_id, "%s:B" % g_id),
                              ("%s:E" % g_id, "%s:E" % f_id)])
            else:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         <-------------
                            g.E           g.B
                """
                if g_l - g_e < 3 or f_l - f_e < 3:
                    continue
                nodes.update(["%s:B" % f_id,
                              "%s:E" % g_id,
                              "%s:B" % g_id,
                              "%s:E" % f_id])
                edges.update([("%s:B" % f_id, "%s:E" % g_id),
                              ("%s:B" % g_id, "%s:E" % f_id)])

    return ig.Graph.DictList(edges=(dict(source=s, target=t) for s, t in edges),
                             vertices=None,
                             directed=True)


def transitive_reduction(sg):
    pass


def draw_graph(sg):
    E = [e.tuple for e in sg.es]
    N = sg.vcount()
    pos = sg.layout('kk')
    
    edge_trace = go.Scatter(x=[i for l in [(pos[s][0], pos[t][0], None) for s, t in E] for i in l],
                            y=[i for l in [(pos[s][1], pos[t][1], None) for s, t in E] for i in l],
                            line=dict(width=0.5, color='black'),
                            mode='lines')
    
    shapes = [make_line(pos[s][0] + (pos[t][0] - pos[s][0]) * 0.7,
                        pos[s][1] + (pos[t][1] - pos[s][1]) * 0.7,
                        pos[t][0],
                        pos[t][1],
                        "black",
                        4,
                        "below")
              for s, t in E]
    
    node_trace = go.Scatter(x=[pos[node][0] for node in range(N)],
                            y=[pos[node][1] for node in range(N)],
                            text=[f"{node['name']}<br>{node.outdegree()} out-nodes" for node in sg.vs],
                            mode='markers',
                            hoverinfo='text',
                            marker=dict(
                                showscale=False,
                                colorscale='YlGnBu',
                                reversescale=True,
                                color=[node.outdegree() for node in sg.vs],
                                size=10,
                                line=dict(width=2)))
    
    layout = go.Layout(width=1000, height=1000,
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       shapes=shapes,
                       hovermode='closest',
                       margin=go.layout.Margin(l=0, r=0, b=0, t=0),
                       showlegend=False)
    fig = go.Figure(data=[edge_trace, node_trace], layout=layout)
    py.iplot(fig)


def construct_unit_graph(overlaps):
    pass
