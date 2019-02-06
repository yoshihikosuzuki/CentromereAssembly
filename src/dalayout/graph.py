from logzero import logger
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt

plt.style.use('ggplot')


def _svs_read_alignment(varmat_i, varmat_j, plot=False):
    n_units_i, n_units_j = varmat_i.shape[0], varmat_j.shape[0]

    def calc_dist_mat(varmat_i, varmat_j):   # TODO: replace to a more efficient code
        dist_mat = np.zeros((n_units_i, n_units_j))
        for i in range(n_units_i):
            for j in range(n_units_j):
                if varmat_i.iloc[i][0] != varmat_j.iloc[j][0]:   # peak, repr, and/or strand is different
                    sim = 0   # i.e. mismatch score for different types of the units = -0.75
                else:
                    assert varmat_i.iloc[i][1].shape[0] == varmat_j.iloc[j][1].shape[0]
                    n_vars = varmat_i.iloc[i][1].shape[0]
                    sim = 1. - float(np.count_nonzero(varmat_i.iloc[i][1] != varmat_j.iloc[j][1])) / n_vars
                dist_mat[i][j] = sim   # similarity
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
    if plot:
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


def svs_read_alignment(read_df_i,
                       read_df_j,
                       varvec_colname,
                       th_n_shared_units,
                       strand,
                       plot=False,
                       debug=False):

    def df_to_composition(df, s):
        return Counter(df.apply(lambda df: (df["peak_id"],
                                            df["repr_id"],
                                            df["strand"] if s == 'f' else 1 - df["strand"]),
                                axis=1))

    def df_to_varmat(df, s):
       return (df if s == 'f'
               else df[::-1]).apply(lambda df: ((df["peak_id"],
                                                 df["repr_id"],
                                                 df["strand"] if s == 'f' else 1 - df["strand"]),
                                                df[varvec_colname]),
                                    axis=1) \
                                    .reset_index(drop=True)

    assert strand in set(['f', 'r']), "Invalid strand value"

    # Check if overlap between the two reads exists at the resolution of representative units
    comp_i = df_to_composition(read_df_i, 'f')
    comp_j = df_to_composition(read_df_j, strand)
    # At least <th_n_shared_units> units must be shared between the two reads
    if sum((comp_i & comp_j).values()) < th_n_shared_units:
        if debug:
            logger.info(f"Composition does not overlap: {comp_i} vs {comp_j}")
        return None

    # Take alignment at the resolution of variant vectors
    varmat_i = df_to_varmat(read_df_i, 'f')   # pd.Series([((peak_id, repr_id, strand), var_vec), ...])
    varmat_j = df_to_varmat(read_df_j, strand)
    ret = _svs_read_alignment(varmat_i, varmat_j, plot=plot)
    if strand == 'r':
        ret["overlap"][3], ret["overlap"][4] = ret["overlap"][5] - ret["overlap"][4], ret["overlap"][5] - ret["overlap"][3]
    return ret


def ava_read_alignment(encodings,
                       reads,
                       varvec_colname="var_vec_global0.0",
                       th_n_shared_units=5,
                       th_align_score=0.5,
                       plot=False):   # TODO: change to a class Overlap?

    def call_svs(strand):
        ret = svs_read_alignment(read_df_i,
                                 read_df_j,
                                 varvec_colname,
                                 th_n_shared_units,
                                 strand)
        if ret is not None and ret["score"] >= th_align_score:
            overlaps.append([read_i, read_j, strand] + ret["overlap"])
        
    read_ids = sorted(set(encodings["read_id"]))
    overlaps = []
    for i, read_i in enumerate(read_ids):
        for read_j in read_ids[i + 1:]:
            # TODO: using only complete units is OK?
            read_df_i = encodings.pipe(lambda df: df[df["read_id"] == read_i]) \
                                 .pipe(lambda df: df[df["type"] == "complete"])
            read_df_j = encodings.pipe(lambda df: df[df["read_id"] == read_j]) \
                                 .pipe(lambda df: df[df["type"] == "complete"])
            if read_df_i.shape[0] == 0 or read_df_j.shape[0] == 0:
                continue
            call_svs('f')
            call_svs('r')
    return overlaps


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
