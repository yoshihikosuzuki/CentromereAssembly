from typing import List
from dataclasses import dataclass, field, InitVar
from multiprocessing import Pool
from logzero import logger
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt
from BITS.utils import run_command, sge_nize, save_pickle

plt.style.use('ggplot')


def _svs_read_alignment(read_i,
                        read_j,
                        strand,
                        varmats,
                        match_th=0.8,
                        indel_penalty=0.3,
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

        if plot:
            # DP matrix only on the optimal alignment path
            dp_opt_path = np.zeros_like(dp)
            dm_opt_path = np.zeros_like(dp)
            for a in alignment:
                dp_opt_path[a[0]][a[1]] = dp[a[0]][a[1]]
                dm_opt_path[a[0]][a[1]] = dist_mat[a[0] - 1][a[1] - 1]
            return (dp, alignment, dp_opt_path, dm_opt_path)
        else:
            return (dp, alignment)

    varmat_i, varmat_j = varmats[(read_i, 'f')], varmats[(read_j, strand)]
    dist_mat = calc_dist_mat()
    if plot:
        dp, alignment, dp_opt_path, dm_opt_path = calc_alignment()
    else:
        dp, alignment = calc_alignment()
    start_i, start_j = alignment[0]
    end_i, end_j = alignment[-1]
    score = dp[end_i][end_j]

    if plot:
        fig = plt.figure(figsize=(25, 10))
        ax1 = fig.add_subplot(141)
        im1 = ax1.imshow(dist_mat, cmap="GnBu", vmin=0.5, vmax=1)
        fig.colorbar(im1)
        ax2 = fig.add_subplot(142)
        im2 = ax2.imshow(dp, cmap="GnBu", vmin=0, vmax=1)
        fig.colorbar(im2)
        ax3 = fig.add_subplot(143)
        im3 = ax3.imshow(dp_opt_path, cmap="GnBu", vmin=0, vmax=1)
        fig.colorbar(im3)
        ax4 = fig.add_subplot(144)
        im4 = ax4.imshow(dm_opt_path, cmap="GnBu", vmin=0.5, vmax=1)
        fig.colorbar(im4)
        fig.show()

    return [start_i, end_i, varmat_i.shape[0], start_j, end_j, varmat_j.shape[0], score, alignment]


def svs_read_alignment(read_i, read_j, strand, comps, varmats, th_n_shared_units, th_align_score):
    if sum((comps[(read_i, 'f')] & comps[(read_j, strand)]).values()) >= th_n_shared_units:
        overlap = _svs_read_alignment(read_i, read_j, strand, varmats)
        if overlap is not None and overlap[6] >= th_align_score:
            if strand == 'r':
                # NOTE: do not expand these assignments
                overlap[3], overlap[4] = overlap[5] - overlap[4], overlap[5] - overlap[3]
            return [read_i, read_j, strand] + overlap
    return None


def svs_read_alignment_mult(list_pairs, comps, varmats, th_n_shared_units, th_align_score):
    return [svs_read_alignment(read_i,
                               read_j,
                               strand,
                               comps,
                               varmats,
                               th_n_shared_units,
                               th_align_score)
            for read_i, read_j, strand in list_pairs]


@dataclass(repr=False, eq=False)
class Overlap:
    encodings: pd.DataFrame
    varvec_colname: str = "var_vec_global0.0"
    th_n_shared_units: int = 5
    th_align_score: float = 0.5

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
                           self.th_align_score)
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
                                                        "alignment")) \
                                    .sort_values(by="strand") \
                                    .sort_values(by="read_j", kind="mergesort") \
                                    .sort_values(by="read_i", kind="mergesort") \
                                    .reset_index(drop=True)


def construct_string_graph(overlaps, th_score=0.5):
    sg = nx.DiGraph()
    for i, overlap in overlaps.iterrows():
        f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l, score, alignment = overlap
        if score < th_score:
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
                sg.add_edge("%s:B" % f_id, "%s:E" % g_id, label=(g_id, g_b, g_l),
                            length=abs(g_b - g_l))
                sg.add_edge("%s:B" % g_id, "%s:E" % f_id, label=(f_id, f_e, f_l),
                            length=abs(f_e - f_l),)

    return sg


def draw_graph(sg, node_size=50, figsize=(20, 20)):
    plt.figure(figsize=figsize)
    pos = nx.spring_layout(sg)
    nx.draw(sg, pos, node_size=node_size)
    nx.draw_networkx_labels(sg, pos)
    plt.show()


def construct_unit_graph(overlaps):
    pass
