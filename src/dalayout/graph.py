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


def _svs_read_alignment(read_i,
                        read_j,
                        strand,
                        varmats,
                        match_th=0.75,
                        indel_penalty=0.2,
                        plot=False):

    def calc_dist_mat():
        return np.array([[0 if varmat_i.iloc[i][0] != varmat_j.iloc[j][0]   # different unit type
                          else 1. - float(np.count_nonzero(varmat_i.iloc[i][1] != varmat_j.iloc[j][1])) / varmat_i.iloc[i][1].shape[0]   # similarity
                          for j in range(varmat_j.shape[0])]
                         for i in range(varmat_i.shape[0])])

    def calc_alignment():
        # Smith-Waterman local alignment
        dp = np.zeros((dist_mat.shape[0] + 1, dist_mat.shape[1] + 1))
        for i in range(1, dp.shape[0]):
            for j in range(1, dp.shape[1]):
                dp[i][j] = max([dp[i - 1][j - 1] + dist_mat[i - 1][j - 1] - match_th,
                                dp[i - 1][j] - indel_penalty,
                                dp[i][j - 1] - indel_penalty])   # TODO: reconsider the scoring system

        # Find the starting point of traceback
        if np.max(dp[-1][1:]) >= np.max(dp.T[-1][1:]):   # maximum is on the end row
            argmax = np.array([dp.shape[0] - 1, np.argmax(dp[-1][1:]) + 1])
        else:   # on the end column
            argmax = np.array([np.argmax(dp.T[-1][1:]) + 1, dp.shape[1] - 1])

        # Traceback
        alignment = [argmax]   # [(pos_i, pos_j), ...] from the traceback starting point
        while True:
            if argmax[0] == 1 or argmax[1] == 1:   # reached the start row or colum
                break
            diag = dp[argmax[0] - 1][argmax[1] - 1] + dist_mat[argmax[0] - 1][argmax[1] - 1] - match_th
            horizontal = dp[argmax[0]][argmax[1] - 1] - indel_penalty
            vertical = dp[argmax[0] - 1][argmax[1]] - indel_penalty
            maximum = np.argmax([diag, horizontal, vertical])   # pointer to the next cell
            if maximum == 0:   # diagonal
                argmax = argmax - 1
            elif maximum == 1:   # horizonal
                argmax = argmax - [0, 1]   # copy object
            else:   # vertical
                argmax = argmax - [1, 0]
            alignment = [argmax] + alignment   # add the next cell at the FRONT of the list

        return (dp, alignment)

    def plot_alignment_mat():
        # NOTE: default heatmap is [x (on y-axis) * y (on x-axis)] for a distance matrix [x * y],
        #       so transpose the matrix.
        
        # Distance matrix
        trace1 = go.Heatmap(z=dist_mat.T,
                           text=np.array([[f"{varmat_i.iloc[ri][0]} vs {varmat_j.iloc[ci][0]}<br>%sim={c:.2}"
                                           for ci, c in enumerate(r)]
                                          for ri, r in enumerate(dist_mat)]).T,
                           hoverinfo="text",
                           colorscale='YlGnBu',
                           #colorscale='Greys',
                           zmin=0.6,
                           zmax=1,
                           reversescale=True,
                           showscale=False)
    
        trace2 = go.Scatter(x=[x[0] - 1 for x in alignment],
                            y=[x[1] - 1 for x in alignment],
                            text=[f"{varmat_i.iloc[x[0] - 1][0]} vs {varmat_j.iloc[x[1] - 1][0]}<br>%sim={dist_mat[x[0] - 1][x[1] - 1]:.2}"
                                  for x in alignment],
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
                           text=np.array([[f"{varmat_i.iloc[ri - 1][0]} vs {varmat_j.iloc[ci - 1][0]}<br>%sim={c:.2}"
                                           if ri * ci != 0 else "0"
                                           for ci, c in enumerate(r)]
                                          for ri, r in enumerate(dp)]).T,
                           hoverinfo="text",
                           colorscale='YlGnBu',
                           #colorscale='Greys',
                           reversescale=True,
                           showscale=False)
    
        trace4 = go.Scatter(x=[x[0] for x in alignment],
                            y=[x[1] for x in alignment],
                            text=[f"{varmat_i.iloc[x[0] - 1][0]} vs {varmat_j.iloc[x[1] - 1][0]}<br>%sim={dp[x[0]][x[1]]:.2}"
                                  for x in alignment],
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

    varmat_i, varmat_j = varmats[(read_i, 'f')], varmats[(read_j, strand)]
    dist_mat = calc_dist_mat()
    dp, alignment = calc_alignment()
    start_i, start_j = alignment[0]
    end_i, end_j = alignment[-1]
    score = dp[end_i][end_j]
    mean_score = float(score) / len(alignment)

    if plot:
        plot_alignment_mat()

    return [start_i, end_i, varmat_i.shape[0], start_j, end_j, varmat_j.shape[0], score, mean_score, alignment]


def svs_read_alignment(read_i, read_j, strand, comps, varmats, th_n_shared_units, th_mean_score):
    if sum((comps[(read_i, 'f')] & comps[(read_j, strand)]).values()) >= th_n_shared_units:
        overlap = _svs_read_alignment(read_i, read_j, strand, varmats)
        if overlap is not None and overlap[7] >= th_mean_score:
            if strand == 'r':
                # NOTE: do not expand these assignments
                overlap[3], overlap[4] = overlap[5] - overlap[4], overlap[5] - overlap[3]
            return [read_i, read_j, strand] + overlap
    return None


def svs_read_alignment_mult(list_pairs, comps, varmats, th_n_shared_units, th_mean_score):
    return [svs_read_alignment(read_i,
                               read_j,
                               strand,
                               comps,
                               varmats,
                               th_n_shared_units,
                               th_mean_score)
            for read_i, read_j, strand in list_pairs]


@dataclass(repr=False, eq=False)
class Overlap:
    encodings: pd.DataFrame
    varvec_colname: str = "var_vec_global0.0"
    th_n_shared_units: int = 10   # TODO: XXX: increase the value as 30 or so to save the computation time !!!
    th_mean_score: float = 0.0

    read_ids: List[int] = field(init=False)
    read_dfs: dict = field(init=False)
    comps: dict = field(init=False)
    varmats: dict = field(init=False)

    def __post_init__(self):
        self.read_dfs = {index: df[df["type"] == "complete"]
                         for index, df
                         in self.encodings[self.encodings["type"] == "complete"].groupby("read_id")}
        self.read_ids = sorted(self.read_dfs.keys())
        self.comps = {(read_id, strand): self.df_to_composition(self.read_dfs[read_id], strand)
                      for read_id in self.read_ids
                      for strand in ['f', 'r']}
        self.varmats = {(read_id, strand): self.df_to_varmat(self.read_dfs[read_id], strand)
                        for read_id in self.read_ids
                        for strand in ['f', 'r']}

    def df_to_composition(self, df, strand):
        return Counter(df.apply(lambda df: (df["peak_id"],
                                            df["repr_id"],
                                            df["strand"] if strand == 'f' else 1 - df["strand"]),
                                axis=1))

    def df_to_varmat(self, df, strand):
        return (df if strand == 'f'
                else df[::-1]).apply(lambda df: ((df["peak_id"],
                                                  df["repr_id"],
                                                  df["strand"] if strand == 'f' else 1 - df["strand"]),
                                                 df[self.varvec_colname]),
                                     axis=1).reset_index(drop=True)

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
        # Finalize
        

    def _ava_read_alignment(self, list_pairs, n_core):
        n_pairs = len(list_pairs)
        n_sub = -(-n_pairs // n_core)
        list_pairs_sub = [(list_pairs[i * n_sub:(i + 1) * n_sub - 1],
                           self.comps,
                           self.varmats,
                           self.th_n_shared_units,
                           self.th_mean_score)
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
                                                        "i_start",
                                                        "i_end",
                                                        "i_len",
                                                        "j_start",
                                                        "j_end",
                                                        "j_len",
                                                        "score",
                                                        "mean_score",
                                                        "alignment")) \
                                    .sort_values(by="strand") \
                                    .sort_values(by="read_j", kind="mergesort") \
                                    .sort_values(by="read_i", kind="mergesort") \
                                    .reset_index(drop=True)


def construct_string_graph(overlaps, th_mean_score=0.0, th_overlap_len=10):
    sg = ig.Graph(directed=True)
    for i, overlap in overlaps.iterrows():
        f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l, score, mean_score, alignment = overlap
        if mean_score < th_mean_score or len(alignment) < th_overlap_len:
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
                if f_b < 3 or g_l - g_e < 3:
                    #print("contained 1")
                    continue
                sg.add_vertices(["%s:B" % g_id, "%s:B" % f_id, "%s:E" % f_id, "%s:E" % g_id])
                sg.add_edge("%s:B" % g_id, "%s:B" % f_id, label=(f_id, f_b, 0),
                            length=abs(f_b - 0))
                sg.add_edge("%s:E" % f_id, "%s:E" % g_id, label=(g_id, g_e, g_l),
                            length=abs(g_e - g_l))
            else:
                """
                     f.B         f.E
                  f  ----------->
                  g         <-------------
                            g.E           g.B
                """
                if f_b < 3 or g_e < 3:
                    #print("contained 2")
                    continue
                sg.add_vertices(["%s:E" % g_id, "%s:B" % f_id, "%s:E" % f_id, "%s:B" % g_id])
                sg.add_edge("%s:E" % g_id, "%s:B" % f_id, label=(f_id, f_b, 0),
                            length=abs(f_b - 0))
                sg.add_edge("%s:E" % f_id, "%s:B" % g_id, label=(g_id, g_e, 0),
                            length=abs(g_e - 0))
        else:
            if g_b < g_e:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         ------------->
                            g.B           g.E
                """
                if g_b < 3 or f_l - f_e < 3:
                    #print("contained 3")
                    continue
                sg.add_vertices(["%s:B" % f_id, "%s:B" % g_id, "%s:E" % g_id, "%s:E" % f_id])
                sg.add_edge("%s:B" % f_id, "%s:B" % g_id, label=(g_id, g_b, 0),
                            length=abs(g_b - 0))
                sg.add_edge("%s:E" % g_id, "%s:E" % f_id, label=(f_id, f_e, f_l),
                            length=abs(f_e - f_l))
            else:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         <-------------
                            g.E           g.B
                """
                if g_l - g_e < 3 or f_l - f_e < 3:
                    #print("contained 4")
                    continue
                sg.add_vertices(["%s:B" % f_id, "%s:E" % g_id, "%s:B" % g_id, "%s:E" % f_id])
                sg.add_edge("%s:B" % f_id, "%s:E" % g_id, label=(g_id, g_b, g_l),
                            length=abs(g_b - g_l))
                sg.add_edge("%s:B" % g_id, "%s:E" % f_id, label=(f_id, f_e, f_l),
                            length=abs(f_e - f_l),)

    return sg


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
