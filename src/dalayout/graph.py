from typing import List
from dataclasses import dataclass, field, InitVar
from logzero import logger
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt

plt.style.use('ggplot')


@dataclass(repr=False, eq=False)
class Overlap:
    encodings: pd.DataFrame
    varvec_colname: str = "var_vec_global0.0"
    th_n_shared_units: int = 5
    th_align_score: float = 0.5
    debug: bool = False

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
                      for strand in ['f', 'r']
                      for read_id in self.read_ids}
        self.varmats = {(read_id, strand): self.df_to_varmat(self.read_dfs[read_id], strand)
                        for strand in ['f', 'r']
                        for read_id in self.read_ids}

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
                                    axis=1) \
                                    .reset_index(drop=True)

    def ava_read_alignment(self, plot=False):
        self.overlaps = {}
        self.index = 0
        for i, read_i in enumerate(self.read_ids):
            for read_j in self.read_ids[i + 1:]:
                self.svs_read_alignment(read_i, read_j, 'f')
                self.svs_read_alignment(read_i, read_j, 'r')

        self.overlaps = pd.DataFrame.from_dict(self.overlaps,
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

    def svs_read_alignment(self, read_i, read_j, strand):
        if sum((self.comps[(read_i, 'f')] & self.comps[(read_j, strand)]).values()) >= self.th_n_shared_units:
            overlap = self._svs_read_alignment(read_i, read_j, strand)
            if overlap is not None and overlap[6] >= self.th_align_score:
                if strand == 'r':
                    # NOTE: do not expand these assignments
                    overlap[3], overlap[4] = overlap[5] - overlap[4], overlap[5] - overlap[3]
                self.overlaps[self.index] = [read_i, read_j, strand] + overlap
                self.index += 1
        else:
            if self.debug:
                logger.info(f"Composition mismatch: [{read_i}] {self.comps[(read_i, 'f')]} vs [{read_j}({strand})] {self.comps[(read_j, 'r')]}")

    def _svs_read_alignment(self, read_i, read_j, strand, plot=False):
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
                    dp[i][j] = max([dp[i - 1][j - 1] + dist_mat[i - 1][j - 1] - 0.75,
                                    dp[i - 1][j] - 0.2,
                                    dp[i][j - 1] - 0.2])   # TODO: reconsider the scoring system

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

        varmat_i, varmat_j = self.varmats[(read_i, 'f')], self.varmats[(read_j, strand)]
        dist_mat = calc_dist_mat()
        dp, dp_optimal_path, alignment = calc_alignment()
        start_i, start_j = alignment[0]
        end_i, end_j = alignment[-1]
        score = dp[end_i][end_j]

        if plot:
            fig = plt.figure(figsize=(18, 10))
            ax1 = fig.add_subplot(131)
            im1 = ax1.imshow(dist_mat, cmap="GnBu", vmin=0.5, vmax=1)
            fig.colorbar(im1)
            ax2 = fig.add_subplot(132)
            im2 = ax2.imshow(dp, cmap="GnBu", vmin=0, vmax=1)
            fig.colorbar(im2)
            ax3 = fig.add_subplot(133)
            im3 = ax3.imshow(dp_optimal_path, cmap="GnBu", vmin=0, vmax=1)
            fig.colorbar(im3)
            fig.show()

        return [start_i, end_i, varmat_i.shape[0], start_j, end_j, varmat_j.shape[0], score, alignment]


def construct_string_graph(overlaps):
    sg = nx.DiGraph()
    for overlap in overlaps:
        f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l = overlap
        print(f_id, g_id, strand, f_b, f_e, f_l, g_b, g_e, g_l )

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
