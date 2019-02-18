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


def plot_alignment_mat(read_sig_i, read_sig_j, score_mat, dp, path):
    # NOTE: default heatmap is [x (on y-axis) * y (on x-axis)] for a distance matrix [x * y],
    #       so transpose the matrix.
    
    # Score matrix
    trace1 = go.Heatmap(z=score_mat.T,
                       text=np.array([[f"{ri}: {read_sig_i[ri][:4]} vs {ci}: {read_sig_j[ci][:4]}<br>%sim={c:.2}"
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
                        text=[f"{x[0]}: {read_sig_i[x[0]][:4]} vs {x[1]}: {read_sig_j[x[1]][:4]}<br>%sim={score_mat[x[0]][x[1]]:.2}"
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
                       text=np.array([[f"{ri - 1}: {read_sig_i[ri - 1][:4]} vs {ci - 1}: {read_sig_j[ci - 1][:4]}<br>%sim={c:.2}"
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
                        text=[f"{x[0]}: {read_sig_i[x[0]][:4]} vs {x[1]}: {read_sig_j[x[1]][:4]}<br>%sim={dp[x[0] + 1][x[1] + 1]:.2}"
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


def calc_score_mat(read_sig_i, read_sig_j, match_th):
    return np.array([[0 if (read_sig_i[i][0] != read_sig_j[j][0]
                            or read_sig_i[i][2] != read_sig_j[j][2])
                      else match_th if (read_sig_i[i][3] != "complete"
                                        or read_sig_j[j][3] != "complete")
                      else 0 if read_sig_i[i][1] != read_sig_j[j][1]
                      else 1. - float(np.count_nonzero(read_sig_i[i][4] != read_sig_j[j][4])) / read_sig_i[i][4].shape[0]   # TODO: count only 1-1 matches
                      for j in range(len(read_sig_j))]
                     for i in range(len(read_sig_i))])


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


def _svs_read_alignment(read_sig_i,
                        read_sig_j,
                        match_th=0.7,
                        indel_penalty=0.2,
                        plot=False):
    """
    Perform encoded read vs encoded read alignment.
    """

    # Calculate match score for every unit pair between read_i and read_j
    score_mat = calc_score_mat(read_sig_i, read_sig_j, match_th)

    # Take alignment by solving a DP which intermediates between NW and SW
    dp, path, score = calc_alignment(score_mat, match_th, indel_penalty)
    mean_score = score / len(path)

    if plot:
        plot_alignment_mat(read_sig_i, read_sig_j, score_mat, dp, path)

    return (mean_score, path)


def svs_read_alignment(read_i,
                       read_j,
                       strand,
                       read_sig_i,
                       read_sig_j,
                       reads,
                       th_mean_score,
                       th_ovlp_len):

    mean_score, path = _svs_read_alignment(read_sig_i, read_sig_j)
    if mean_score < th_mean_score:   # Be careful that score depends on <match_th> and <indel_penalty>
        return None

    #logger.debug(f"{read_i} vs {read_j}({strand})")
    overlap_len = sum([max(read_sig_i[i][7], read_sig_j[j][7]) for i, j in path])   # in bp
    if overlap_len < th_ovlp_len:
        return None
    # TODO: overlap length not by sum of the unit lengths but end_bp - start_bp
    # TODO: check consistency of end_bp - start-bp between read i and j

    n_unit_i, n_unit_j = len(read_sig_i), len(read_sig_j)
    start_unit_i, start_unit_j = path[0]
    end_unit_i, end_unit_j = path[-1]

    start_bp_i, end_bp_i = read_sig_i[start_unit_i][5], read_sig_i[end_unit_i][6]
    start_bp_j = read_sig_j[start_unit_j if strand == 0 else end_unit_j][5]
    end_bp_j = read_sig_j[end_unit_j if strand == 0 else start_unit_j][6]
    
    if strand == 1:
        start_unit_j, end_unit_j = n_unit_j - 1 - end_unit_j, n_unit_j - 1 - start_unit_j

    read_len_i, read_len_j = reads.loc[read_i, "length"], reads.loc[read_j, "length"]
    if start_unit_i > 3:
        if strand == 0:
            """
                 f.B         f.E
              f  ----------->
              g         ------------->
                        g.B           g.E
            """
            if n_unit_j - end_unit_j - 1 < 3:
                overlap_type = "contains"
            else:
                overlap_type = "suffix-prefix"
        else:
            """
                 f.B         f.E
              f  ----------->
              g         <-------------
                        g.E           g.B
            """
            if start_unit_j < 3:
                overlap_type = "contains"
            else:
                overlap_type = "suffix-suffix"
    else:
        if n_unit_i - end_unit_i - 1 < 3:
            overlap_type = "contained"
        else:
            if strand == 0:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         ------------->
                            g.B           g.E
                """
                if start_unit_j < 3:
                    overlap_type = "contains"
                else:
                    overlap_type = "prefix-suffix"
            else:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         <-------------
                            g.E           g.B
                """
                if n_unit_j - end_unit_j - 1 < 3:
                    overlap_type = "contains"
                else:
                    overlap_type = "prefix-prefix"
    # TODO: check consistency of overlap type in units and bps between read i and read j
    # (i.e. is the overlap really dovetail even with bp coordinate system?)

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
            read_len_i,
            start_bp_j,
            end_bp_j,
            read_len_j,
            overlap_len,
            overlap_type,
            mean_score,
            path]


def svs_read_alignment_mult(list_pairs,
                            read_sigs,
                            read_comps,
                            th_n_shared_units,
                            th_mean_score,
                            th_ovlp_len):

    return [svs_read_alignment(read_i,
                               read_j,
                               strand,
                               read_sigs[(read_i, 0)],
                               read_sigs[(read_j, strand)],
                               th_mean_score,
                               th_ovlp_len)
            for read_i, read_j, strand in list_pairs
            if (sum((read_comps[(read_i, 0)]
                     & read_comps[(read_j, strand)]).values())
                >= th_n_shared_units)]   # filter by representative units composition


@dataclass(repr=False, eq=False)
class Overlap:
    encodings: InitVar[pd.DataFrame]
    varvec_colname: str = "var_vec_global0.15"
    th_n_shared_units: int = 10   # TODO: change to bp
    th_mean_score: float = 0.04
    th_ovlp_len: int = 3000   # in bp

    read_sigs: dict = field(init=False)
    read_comps: dict = field(init=False)

    def __post_init__(self, encodings):
        gb = encodings.groupby("read_id")

        # To quickly iterate over the units in a read during alignment, convert necessary data in encodings
        # into a list for each read and for each strand
        self.read_sigs = {(read_id, 0): list(df.apply(lambda df: (df["peak_id"],
                                                                  df["repr_id"],
                                                                  df["strand"],
                                                                  df["type"],
                                                                  df[self.varvec_colname],
                                                                  df["start"],
                                                                  df["end"],
                                                                  df["length"]),
                                                      axis=1))
                          for read_id, df in gb}
        self.read_sigs.update({(k[0], 1): [(x[0],
                                            x[1],
                                            (x[2] + 1) % 2,   # revcomp
                                            x[3],
                                            x[4],
                                            x[5],
                                            x[6],
                                            x[7])   # TODO: can we convert start/end positions here?
                                           for x in reversed(v)]
                               for k, v in self.read_sigs.items()})

        # Composition of representative units for each read and for each strand
        # Used for initial filtering of alignment candidate
        self.read_comps = {k: Counter([x[:3] for x in v])
                           for k, v in self.read_sigs.items()}

    def svs_read_alignment(self, read_i, read_j, strand, match_th=0.7, indel_penalty=0.2, plot=False):
        return _svs_read_alignment(self.read_sigs[(read_i, 0)],
                                   self.read_sigs[(read_j, strand)],
                                   match_th,
                                   indel_penalty,
                                   plot)

    def list_up_pairs(self):
        """
        List up pairs of reads whose alignment will be taken.
        Only pairs are filtered by the composition of representative units in reads.
        """

        read_ids = sorted({x[0] for x in self.read_comps.keys()})
        return [(read_i, read_j, strand)
                for i, read_i in enumerate(read_ids)
                for read_j in read_ids[i + 1:]
                for strand in [0, 1]]

    def ava_read_alignment(self, n_core):
        """
        Non-distributed version of all-vs-all read overlap calculation.
        """

        self._ava_read_alignment(self.list_up_pairs(), n_core)

    def ava_read_alignment_distribute(self, n_distribute, n_core):
        """
        Distributed all-vs-all read overlap calculation.
        """

        # Split the tasks
        list_pairs = self.list_up_pairs()
        n_sub = -(-len(list_pairs) // n_distribute)
        list_pairs_sub = [list_pairs[i * n_sub:(i + 1) * n_sub - 1]
                          for i in range(n_distribute)]

        # Prepare data for each distributed job and submit them
        save_pickle(self, "overlap_obj.pkl")
        n_digit = int(np.log10(n_distribute) + 1)
        for i, lp in enumerate(list_pairs_sub):
            index = str(i + 1).zfill(n_digit)
            save_pickle(lp, f"list_pairs.{index}.pkl")
            script_fname = f"ava_pair.sge.{index}"
            with open(script_fname, 'w') as f:
                f.write(sge_nize(' '.join([f"run_ava_sub.py",
                                           f"-n {n_core}",
                                           f"overlap_obj.pkl",
                                           f"list_pairs.{index}.pkl",
                                           f"overlaps.{index}.pkl"]),
                                 job_name="run_ava_sub",
                                 n_core=n_core,
                                 wait=False))
            run_command(f"qsub {script_fname}")

    def _ava_read_alignment(self, list_pairs, n_core):
        n_sub = -(-len(list_pairs) // n_core)
        list_pairs_sub = [(list_pairs[i * n_sub:(i + 1) * n_sub - 1],
                           self.read_sigs,
                           self.read_comps,
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
                                                        "overlap_len",
                                                        "mean_score",
                                                        "path")) \
                                    .sort_values(by="strand") \
                                    .sort_values(by="read_j", kind="mergesort") \
                                    .sort_values(by="read_i", kind="mergesort") \
                                    .reset_index(drop=True)


def construct_string_graph(overlaps, reads, th_mean_score=0.04, th_overlap_len=3000):
    # NOTE: Because igraph prefers static graph construction, first list the vertices and edges up.
    nodes, edges = set(), set()
    for i, overlap in overlaps.iterrows():
        f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l, f_bb, f_be, f_bl, g_bb, g_be, g_bl, ovlp_len, ovlp_type, mean_score, path = overlap

        if mean_score < th_mean_score or ovlp_len < th_overlap_len:
            continue

        if ovlp_type in ["contains", "contained"]:
            continue
        elif ovlp_type == "suffix-prefix":
            nodes.update(["%s:B" % g_id,
                          "%s:B" % f_id,
                          "%s:E" % f_id,
                          "%s:E" % g_id])
            edges.update([("%s:B" % g_id, "%s:B" % f_id, f_bb),
                          ("%s:E" % f_id, "%s:E" % g_id, g_bl - g_be)])
        elif ovlp_type == "suffix-suffix":
            nodes.update(["%s:E" % g_id,
                          "%s:B" % f_id,
                          "%s:E" % f_id,
                          "%s:B" % g_id])
            edges.update([("%s:E" % g_id, "%s:B" % f_id, f_bb),
                          ("%s:E" % f_id, "%s:B" % g_id, g_bb)])
        elif ovlp_type == "prefix-suffix":
            nodes.update(["%s:B" % f_id,
                          "%s:B" % g_id,
                          "%s:E" % g_id,
                          "%s:E" % f_id])
            edges.update([("%s:B" % f_id, "%s:B" % g_id, g_bb),
                          ("%s:E" % g_id, "%s:E" % f_id, f_bl - f_be)])
        else:   # prefix-prefix
            nodes.update(["%s:B" % f_id,
                          "%s:E" % g_id,
                          "%s:B" % g_id,
                          "%s:E" % f_id])
            edges.update([("%s:B" % f_id, "%s:E" % g_id, g_bl - g_be),
                          ("%s:B" % g_id, "%s:E" % f_id, f_bl - f_be)])

    return ig.Graph.DictList(edges=(dict(source=s, target=t, length=l) for s, t, l in edges),
                             vertices=None,
                             directed=True)


def transitive_reduction(sg):
    v_mark = ["vacant" for v in sg.vs]
    e_reduce = {e.tuple: False for e in sg.es}
    FUZZ = 10   # TODO: in bp; assuming no unit shifts

    for v in sg.vs:
        if v.outdegree() == 0:
            continue

        oes = sorted(sg.es.select(_source=v.index), key=lambda x: x["length"])
        longest = oes[-1]["length"] + FUZZ
        for oe in oes:
            v_mark[oe.target] = "inplay"

        for oe in oes:
            if v_mark[oe.target] == "inplay":
                ooes = sorted(sg.es.select(_source=oe.target), key=lambda x: x["length"])
                for ooe in ooes:
                    if oe["length"] + ooe["length"] <= longest and v_mark[ooe.target] == "inplay":
                        v_mark[ooe.target] = "eliminated"

        for oe in oes:
            ooes = sorted(sg.es.select(_source=oe.target), key=lambda x: x["length"])
            if len(ooes) > 1:
                shortest = ooes[0].target
                if v_mark[shortest] == "inplay":
                    v_mark[shortest] == "eliminated"
            for ooe in ooes:
                if ooe["length"] < FUZZ and v_mark[ooe.target] == "inplay":
                    v_mark[ooe.target] = "eliminated"

        for oe in oes:
            if v_mark[oe.target] == "eliminated":
                e_reduce[oe.tuple] = True   # TODO: confirm revcomp edges will be also removed in the same way
            v_mark[oe.target] = "vacant"

    # Re-construct a graph
    return ig.Graph.DictList(edges=(dict(source=e["source"],
                                         target=e["target"],
                                         length=e["length"])
                                    for e in sg.es
                                    if not e_reduce[e.tuple]),
                             vertices=None,
                             directed=True)


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
