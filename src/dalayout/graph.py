from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
import numpy as np
import pandas as pd
import networkx as nx
from interval import interval
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
from BITS.seq import revcomp
from BITS.utils import run_command
from BITS.run import run_edlib
import consed

plt.style.use('ggplot')


def _svs_read_alignment(varmat_i, varmat_j, plot=False):
    n_units_i, n_units_j = varmat_i.shape[0], varmat_j.shape[0]

    def calc_dist_mat(varmat_i, varmat_j):   # TODO: replace to a more efficient code
        dist_mat = np.zeros((n_units_i, n_units_j))
        assert varmat_i.shape[1] == varmat_j.shape[1]
        n_vars = varmat_i.shape[1]
        for i in range(n_units_i):
            for j in range(n_units_j):
                dist_mat[i][j] = 1. - float(np.count_nonzero(varmat_i[i] != varmat_j[j])) / n_vars   # similarity
        return dist_mat

    dist_mat = calc_dist_mat(varmat_i, varmat_j)

    def calc_alignment(dist_mat):
        dp = np.zeros((dist_mat.shape[0] + 1, dist_mat.shape[1] + 1))

        # Smith-Waterman local alignment
        for i in range(1, dp.shape[0]):
            for j in range(1, dp.shape[1]):
                dp[i][j] = max([dp[i - 1][j - 1] + dist_mat[i - 1][j - 1] - 0.75,
                                dp[i - 1][j] - 0.2,
                                dp[i][j - 1] - 0.2])   # TODO: make the scores more sophisticated
        
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
            diag = dp[argmax[0] - 1][argmax[1] - 1] + dist_mat[argmax[0] - 1][argmax[1] - 1] - 0.75
            horizontal = dp[argmax[0]][argmax[1] - 1] - 0.2
            vertical = dp[argmax[0] - 1][argmax[1]] - 0.2
            maximum = np.argmax([diag, horizontal, vertical])   # pointer to the next cell
            if maximum == 0:   # diagonal
                argmax = argmax - 1
            elif maximum == 1:   # horizonal
                argmax = argmax - [0, 1]   # copy object
            else:   # vertical
                argmax = argmax - [1, 0]
            alignment = [argmax] + alignment   # add the next cell at the FRONT of the list
    
        # DP matrix only on the optimal alignment path
        dp_optimal_path = np.zeros_like(dp)
        for a in alignment:
            dp_optimal_path[a[0]][a[1]] = dp[a[0]][a[1]]

        return (dp, dp_optimal_path, alignment)

    dp, dp_optimal_path, alignment = calc_alignment(dist_mat)

    start_i, start_j = alignment[0]
    end_i, end_j = alignment[-1]
    overlap = [start_i, end_i, n_units_i, start_j, end_j, n_units_j]
        
    score = dp[end_i][end_j]   # TODO: change the score
    if plot and score >= 0.5:
        fig = plt.figure(figsize=(18, 10))
        ax1 = fig.add_subplot(131)
        #ax1.set_title(f"{f_id} vs {g_id}({g_s}): dist mat")
        im1 = ax1.imshow(dist_mat, cmap="GnBu", vmin=0.5, vmax=1)
        fig.colorbar(im1)
        ax2 = fig.add_subplot(132)
        #ax2.set_title(f"{f_id} vs {g_id}({g_s}): dp mat")
        im2 = ax2.imshow(dp, cmap="GnBu", vmin=0, vmax=1)
        fig.colorbar(im2)
        ax3 = fig.add_subplot(133)
        #ax3.set_title(f"{f_id} vs {g_id}({g_s}): opt path")
        im3 = ax3.imshow(dp_optimal_path, cmap="GnBu", vmin=0, vmax=1)
        fig.colorbar(im3)
        fig.show()

    return {"alignment": alignment,
            "score": score,
            "overlap": overlap}


def svs_read_alignment(read_df_i, read_df_j, read_len_i, read_len_j, plot_dp=False):
    def df_to_varmat(df):
        return np.array(list(df["var_vec"].values))

    # TODO: XXX: before converting to a variant matrix, check the strand of the units!!!
    # or, keep the strand information (and repr_id information) in addition to variant vector
    read_id_i, read_id_j = read_df_i.iloc[0]["read_id"], read_df_j.iloc[0]["read_id"]
    varmat_i, varmat_j = df_to_varmat(read_df_i), df_to_varmat(read_df_j)
    ret_f = _svs_read_alignment(varmat_i, varmat_j)
    ret_rc = _svs_read_alignment(varmat_i,  varmat_j[::-1])
    ret, strand = (ret_f, "f") if ret_f["score"] >= ret_rc["score"] else (ret_rc, "r")
    if strand == "r":
        ret.overlap[6], ret.overlap[7] = ret.overlap[8] - ret.overlap[7], ret.overlap[8] - ret.overlap[6]
    return [read_id_i, read_id_j, strand] + ret["overlap"]


def ava_read_alignment(encodings, reads):
    return pd.concat([svs_read_alignment(i,
                                         j,
                                         reads.loc[read_i, "length"],
                                         reads.loc[read_j, "length"])
                      for i, j in pairs])


def plot_read_as_varmat(varmat, figsize=(10, 10)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(varmat)
    fig.show()


def construct_string_graph(overlaps):
    sg = nx.DiGraph()
    for overlap in overlaps:
        f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l = overlap
        print(f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l )

        if strand == "r":  # revered alignment, swapping the begin and end coordinates
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
                    print("contained 1")
                    continue
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
                    print("contained 2")
                    continue
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
                    print("contained 3")
                    continue
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
                    print("contained 4")
                    continue
                sg.add_edge("%s:B" % f_id, "%s:E" % g_id, label=(g_id, g_b, g_l),
                            length=abs(g_b - g_l))
                sg.add_edge("%s:B" % g_id, "%s:E" % f_id, label=(f_id, f_e, f_l),
                            length=abs(f_e - f_l),)
    

    # drawing 1
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(sg)
    nx.draw(sg, pos, node_size=500)
    nx.draw_networkx_labels(sg, pos)
    plt.show()

    # drawing 2
    nx.draw_networkx(sg)
    # TODO: another drawing?


def construct_unit_graph(overlaps):
    pass
