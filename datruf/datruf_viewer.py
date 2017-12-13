import os
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib.image as img
import plotly.offline as py
import plotly.graph_objs as go

from datruf_io import (load_tr_intervals,
                       load_alignments,
                       load_paths)

from datruf_core import (calc_cover_set,
                         calc_min_cover_set,
                         take_consensus,
                         trace_path)

from datruf_utils import (run_command,
                          make_line)

plt.style.use('ggplot')
py.init_notebook_mode()


class Viewer:
    """
    Dot plot/alignment plot/alignment path plot viewer for one read.
    The 'read_id' (int) is same as the one used in DAZZ_DB.
    This viewer is assumed to be used in Jupyter Notebook for now.
    """

    def __init__(self, root_dir, db_file, las_file, out_dir, gepard):
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        self.root_dir = root_dir
        self.db_file = db_file
        self.las_file = las_file
        self.out_dir = out_dir
        self.gepard = gepard

    def _calc_cover_set(self):
        if not hasattr(self, 'tr_intervals'):
            # [(start, end), ...]
            self.tr_intervals = load_tr_intervals(self)

        if not hasattr(self, 'alignments'):
            # pd.DataFrame([abpos, aepos, bbpos, bepos(, distance)])
            self.alignments = load_alignments(self)

        if not hasattr(self, 'cover_set'):
            # set((bb, ae, ab, be), ...)
            self.cover_set = calc_cover_set(self)

        if not hasattr(self, 'min_cover_set'):
            # set((ab, ae, bb, be), ...)
            self.min_cover_set = calc_min_cover_set(self.cover_set)

    def _load_paths(self):
        self._calc_cover_set()

        if not hasattr(self, 'paths'):
            # [(ab, ae, bb, be, aseq, bseq, symbols), ...]
            self.paths = load_paths(self)

    def set_read(self, read_id):
        # Initialize
        for attr in ['tr_intervals',
                     'alignments',
                     'cover_set',
                     'min_cover_set',
                     'paths']:
            if hasattr(self, attr):
                delattr(self, attr)

        # Set a new read
        self.read_id = int(read_id)

        # Calculate read length
        command = ("DBdump -s %s %d | awk '$1 == \"S\" {print $2}'"
                   % (self.db_file, self.read_id))
        out = run_command(command)
        self.read_len = int(out.strip())

    def dot_plot(self):
        out_fasta = self.out_dir + str(self.read_id) + ".fasta"
        out_dotplot = self.out_dir + str(self.read_id) + ".png"

        # Generate fasta
        if not os.path.isfile(out_fasta):
            command = ("DBshow %s %d > %s"
                       % (self.db_file, self.read_id, out_fasta))
            run_command(command)

        # Calculate dot plot
        if not os.path.isfile(out_dotplot):
            command = ("unset DISPLAY;"   # no usage of any display
                       "%s -seq1 %s -seq2 %s -outfile %s"
                       % (self.gepard, out_fasta, out_fasta, out_dotplot))
            run_command(command)

        # Show dot plot
        fig, ax = plt.subplots(figsize=(11, 11))
        ax.tick_params(labelbottom="off", bottom="off")
        ax.tick_params(labelleft="off", left="off")
        # this assignment and plt.show() are necessary to show only one figure
        tmp = plt.imshow(img.imread(out_dotplot))
        plt.show()

    def alignment_plot(self, show_grid=False):
        self._calc_cover_set()

        # Put TR intervals reported by datander on diagonal
        # diagonal
        shapes = [make_line(0, 0, self.read_len, self.read_len, 'yellow', 3)]
        for start, end in self.tr_intervals:
            # TR interval
            shapes.append(make_line(start, start, end, end, 'black', 3))

        # Determine color of each alignment
        for alignment in self.alignments.iterrows():
            ab, ae, bb, be = alignment[1][["abpos", "aepos", "bbpos", "bepos"]]

            if show_grid is True:
                # Grid lines from all start/end points of the alignments
                col = 'green'
                width = 0.2
                shapes.append(make_line(ab, 0, ab, self.read_len, col, width))
                shapes.append(make_line(ae, 0, ae, self.read_len, col, width))
                shapes.append(make_line(0, bb, self.read_len, bb, col, width))
                shapes.append(make_line(0, be, self.read_len, be, col, width))

            if (ab, ae, bb, be) in self.min_cover_set:
                shapes.append(make_line(ab, bb, ae, be, 'purple', 3))
                continue

            if (bb, ae, ab, be) in self.cover_set:
                col = 'red'
            elif not (0.95 <= float(ae - ab) / (be - bb) <= 1.05):
                col = 'yellow'   # abnormal slope
            elif ab - be <= 20:
                col = 'blue'   # TR alignment but not in cover set
            else:
                col = 'black'   # others

            shapes.append(make_line(ab, bb, ae, be, col, 1))

        # Show alignment plot
        trace1 = go.Scatter(x=list(self.alignments["abpos"]),
                            y=list(self.alignments["bbpos"]),
                            text=list(range(1, len(self.alignments) + 1)),
                            textposition="top right",
                            textfont=dict(size=10, color="red"),
                            mode='text',
                            name="alignment #")

        trace2 = go.Scatter(x=list(self.alignments["abpos"]),
                            y=list(self.alignments["bbpos"]),
                            mode='markers',
                            name="start")

        trace3 = go.Scatter(x=list(self.alignments["aepos"]),
                            y=list(self.alignments["bepos"]),
                            mode='markers',
                            name="end")

        trace4 = go.Scatter(x=[x[0] for x in self.tr_intervals],
                            y=[x[0] for x in self.tr_intervals],
                            text=list(range(1, len(self.tr_intervals) + 1)),
                            textposition="top right",
                            textfont=dict(size=10, color="black"),
                            mode="text",
                            name="interval #")

        layout = go.Layout(width=725, height=725,
                           hovermode='closest',
                           xaxis=dict(showgrid=False,
                                      range=[0, self.read_len + 100]),
                           yaxis=dict(showgrid=False,
                                      range=[0, self.read_len],
                                      autorange='reversed'),
                           shapes=shapes)

        py.iplot(go.Figure(data=[trace1, trace2, trace3, trace4],
                           layout=layout))

    def path_plot(self, snake=True):   # TODO: make a class for path
        self._load_paths()

        # diagonal
        shapes = [make_line(0, 0, self.read_len, self.read_len, 'yellow', 3)]
        # path trace
        self.divided_paths = []
        for path in self.paths:
            divided_path, add_shapes = trace_path(path, plot=True, snake=snake)
            self.divided_paths.append(divided_path)
            shapes.extend(add_shapes)

        # Show path(s) (and reflecting snake(s))
        trace1 = go.Scatter(x=[x[0] for x in self.paths],
                            y=[x[2] for x in self.paths],
                            mode='markers',
                            name="start")

        trace2 = go.Scatter(x=[x[1] for x in self.paths],
                            y=[x[3] for x in self.paths],
                            mode='markers',
                            name="end")

        layout = go.Layout(width=1000, height=1000,
                           hovermode='closest',
                           xaxis=dict(showgrid=False,
                                      range=[0, self.read_len + 100]),
                           yaxis=dict(showgrid=False,
                                      range=[0, self.read_len],
                                      autorange='reversed'),
                           shapes=shapes)

        py.iplot(go.Figure(data=[trace1, trace2], layout=layout))

    def consensus(self):
        self._load_paths()
        
        if not hasattr(self, 'divided_paths'):
            self.divided_paths = []
            for path in self.paths:
                divided_path, dist = trace_path(path, plot=False)
                self.divided_paths.append(divided_path)

        for path in self.divided_paths:
            DAG = take_consensus(*path[4:])

            # Show DAG
            plt.figure(figsize=(18, 10))
            plt.axis("off")
            #pos = nx.spectral_layout(DAG)
            #pos = nx.circular_layout(DAG)
            #pos = graphviz_layout(DAG, prog="dot")
            pos = graphviz_layout(DAG, prog="neato")
            edge_weights = nx.get_edge_attributes(DAG, 'weight')
            nx.draw_networkx(DAG, pos, with_labels=False, node_size=1, font_size=1)   # TODO: output as dot file
            #nx.draw_networkx_edge_labels(DAG, pos, edge_labels=edge_weights)
