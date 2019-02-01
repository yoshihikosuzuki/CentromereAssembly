from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
import numpy as np
import pandas as pd
from BITS.run import run_edlib
from BITS.utils import print_log, NoDaemonPool
import consed


def encoding_to_read_matrix(read_id, encodings, repr_units):
    fname_prefix = f"peak_{peak_id}_repr_{repr_id}.raw_units"



read_matrix = pd.DataFrame(units).groupby(0)[4].apply(lambda s: np.array(list(s.as_matrix())))


# dist mat heatmap
for mat in read_matrix:
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(mat)
    fig.show()


def dotplot_matrix(m, n):
    ret = np.zeros((m.shape[0], n.shape[0]))
    assert m.shape[1] == n.shape[1]
    L = m.shape[1]
    for i in range(m.shape[0]):
        for j in range(n.shape[0]):
            ret[i][j] = 1. - float(np.count_nonzero(m[i] != n[j])) / L   # similarity
    return ret


def find_alignment(dist_mat, f_id, g_id, g_s):
    dp = np.zeros((dist_mat.shape[0] + 1, dist_mat.shape[1] + 1))
    for i in range(1, dp.shape[0]):
        for j in range(1, dp.shape[1]):
            #print(f"@({i},{j}) {dp[i - 1][j - 1] + dist_mat[i - 1][j - 1] - 0.7}, {dp[i - 1][j] - 0.1}, {dp[i][j - 1] - 0.1}")
            dp[i][j] = max([dp[i - 1][j - 1] + dist_mat[i - 1][j - 1] - 0.75, dp[i - 1][j] - 0.2, dp[i][j - 1] - 0.2])   # TODO: need indel penalty?
    
    # find the starting point of traceback
    if np.max(dp[-1][1:]) >= np.max(dp.T[-1][1:]):   # maximum is on the bottom row
        argmax = np.array([dp.shape[0] - 1, np.argmax(dp[-1][1:]) + 1])
    else:
        argmax = np.array([np.argmax(dp.T[-1][1:]) + 1, dp.shape[1] - 1])
    #print(argmax)
    alignment = [argmax]
    while True:   # tracebacak
        if argmax[0] == 1 or argmax[1] == 1:
            break
        diag = dp[argmax[0] - 1][argmax[1] - 1] + dist_mat[argmax[0] - 1][argmax[1] - 1] - 0.75
        horizontal = dp[argmax[0]][argmax[1] - 1] - 0.2
        vertical = dp[argmax[0] - 1][argmax[1]] - 0.2
        maximum = np.argmax([diag, horizontal, vertical])
        if maximum == 0:
            argmax = argmax - 1
        elif maximum == 1:
            argmax = argmax - [0, 1]   # copy object
        else:
            argmax = argmax - [1, 0]
        alignment = [argmax] + alignment

    #print(alignment)
    # backtraced optimal alignment path
    path = np.zeros_like(dp)
    for a in alignment:
        path[a[0]][a[1]] = dp[a[0]][a[1]]
        
    overlap = [alignment[0][0], alignment[-1][0], dp.shape[0] - 1, alignment[0][1], alignment[-1][1], dp.shape[1] - 1]
        
    score = dp[alignment[-1][0]][alignment[-1][1]]
    if score >= 0.5:
        # show score matrix and dp matrix as heatmap
        #"""
        fig = plt.figure(figsize=(18, 10))
        ax1 = fig.add_subplot(131)
        ax1.set_title(f"{f_id} vs {g_id}({g_s}): dist mat")
        im1 = ax1.imshow(dist_mat, cmap="GnBu", vmin=0.5, vmax=1)
        fig.colorbar(im1)
        ax2 = fig.add_subplot(132)
        ax2.set_title(f"{f_id} vs {g_id}({g_s}): dp mat")
        im2 = ax2.imshow(dp, cmap="GnBu", vmin=0, vmax=1)
        fig.colorbar(im2)
        ax3 = fig.add_subplot(133)
        ax3.set_title(f"{f_id} vs {g_id}({g_s}): opt path")
        im3 = ax3.imshow(path, cmap="GnBu", vmin=0, vmax=1)
        fig.colorbar(im3)
        fig.show()
        #"""
        return (alignment, score, overlap)

    return (alignment, score)


overlaps = []

for i in range(read_matrix.shape[0] - 1):
    scores = []
    scores_rc = []
    for j in range(i + 1, read_matrix.shape[0]):
        ret = find_alignment(dotplot_matrix(read_matrix.iloc[i], read_matrix.iloc[j]), i, j, "f")
        if len(ret) == 2:
            alignment, score = ret
        else:
            alignment, score, overlap = ret
            overlap = [i, j, "f"] + overlap
            overlaps.append(overlap)
        ret = find_alignment(dotplot_matrix(read_matrix.iloc[i], read_matrix.iloc[j][::-1]), i, j, "r")
        if len(ret) == 2:
            alignment_rc, score_rc = ret
        else:
            alignment_rc, score_rc, overlap = ret
            overlap = [i, j, "r"] + overlap
            overlap[6], overlap[7] = overlap[8] - overlap[7], overlap[8] - overlap[6]
            overlaps.append(overlap)
        #print(alignment)
        #print(score)
    #fig, ax = plt.subplots()
    #ax.hist(scores, bins=50)
    #ax.hist(scores_rc, bins=50)
    #plt.show()



import networkx as nx


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


plt.figure(figsize=(12, 12))
pos = nx.spring_layout(sg)
nx.draw(sg, pos, node_size=500)
nx.draw_networkx_labels(sg, pos)
plt.show()


nx.draw_networkx(sg)   # another drawing?



# TODO: write experimental "unit graph" here


