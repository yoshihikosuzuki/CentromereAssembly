import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go
from logzero import logger
from BITS.utils import run_command, make_line
from .load import load_tr_intervals, load_alignments, load_paths
from .core import calc_cover_set, calc_min_cover_set

plt.style.use('ggplot')


class Viewer:
    """
    Dot plot/alignment plot/alignment path plot viewer for one read.
    The 'read_id' (int) is same as the one used in DAZZ_DB.
    This viewer is assumed to be used in Jupyter Notebook for now.

    Basic usage:
        ```
        py.init_notebook_mode(connected=True)
        from datruf import Viewer
        v = Viewer(db_file, las_file, out_dir, gepard_command)
        v.show(read_id, path_plot=True, show_dag=True)
    """

    def __init__(self, db_file, las_file, out_dir, gepard, encodings=None):
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        self.db_file = db_file
        self.las_file = las_file
        self.out_dir = out_dir
        self.gepard = gepard
        self.encodings = encodings

    def _calc_cover_set(self):
        if not hasattr(self, 'tr_intervals'):
            # [(start, end), ...]
            self.tr_intervals = load_tr_intervals(self)

        if not hasattr(self, 'alignments'):
            # pd.DataFrame([abpos, aepos, bbpos, bepos, distance])
            # sorted by distance -> abpos
            self.alignments = load_alignments(self)

        if not hasattr(self, 'cover_set'):
            # set((bb, ae, ab, be), ...)
            self.cover_set = calc_cover_set(self)

        if not hasattr(self, 'min_cover_set'):
            # set((ab, ae, bb, be), ...)
            self.min_cover_set = calc_min_cover_set(self.cover_set)

    def _load_paths(self, plot=False, show_snake=False):
        self._calc_cover_set()   # load_paths needs min_cover_set

        if not hasattr(self, 'paths'):
            # [Path, ...]
            self.paths = load_paths(self)   # only alignments are loaded

        for path in self.paths:
            if plot is True:
                if not hasattr(path, 'shapes'):
                    path.split_alignment(plot=True, snake=show_snake)
            else:
                if not hasattr(path, 'unit_alignments'):
                    path.split_alignment()

    def show(self,
             read_id,
             dot_plot=False,
             varmat_plot=False,
             alignment_plot=True,
             show_grid=False,   # for alignment plot
             path_plot=False,
             show_snake=True,   # for alignment path plot
             save_html=False,   # instead of drawing in the current cell
             show_dag=False):
        """
        All-in-one drawing function for a single read.
        By default, this method generates only dot plot and alignment plot.
        You can additionally plot alignment path plot and DAG of the units.
        """

        self.set_read(read_id)

        if dot_plot:
            self.dot_plot()
        if varmat_plot:
            self.varmat_plot()
        if alignment_plot:
            self.alignment_plot(show_grid, save_html)
        if path_plot:
            self.path_plot(show_snake, save_html)
        if show_dag:
            self.show_dag()

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
        command = (f"DBdump -s {self.db_file} {self.read_id} "
                   f"| awk '$1 == \"S\" {{print $2}}'")
        out = run_command(command)
        self.read_len = int(out.strip())

    def dot_plot(self):
        out_fasta = self.out_dir + str(self.read_id) + ".fasta"
        out_dotplot = self.out_dir + str(self.read_id) + ".png"

        # Generate fasta
        if not os.path.isfile(out_fasta):
            run_command(f"DBshow {self.db_file} {self.read_id} > {out_fasta}")

        # Calculate dot plot
        if not os.path.isfile(out_dotplot):
            command = (f"unset DISPLAY;"   # no usage of any display
                       f"{self.gepard} -seq1 {out_fasta} -seq2 {out_fasta} -outfile {out_dotplot}")
            run_command(command)

        # Show dot plot
        fig, ax = plt.subplots(figsize=(11, 11))
        ax.tick_params(labelbottom=False, bottom=False)
        ax.tick_params(labelleft=False, left=False)
        # this assignment and plt.show() are necessary to show only one figure
        tmp = plt.imshow(img.imread(out_dotplot))
        plt.show()

    def varmat_plot(self, colname="var_vec_global0.0"):   # TODO: change to integrated var. vec. viewer
        if self.encodings is None:
            logger.info("Encodings are not specified. Exit.")
            return
        read_df = self.encodings.pipe(lambda df: df[df["read_id"] == self.read_id]) \
                                .pipe(lambda df: df[df["type"] == "complete"])
        if read_df.shape[0] == 0:
            logger.info("The read does not have any complete variant vector unit. Exit.")
            return
        varmat = np.array(list(read_df[colname].values))
        plt.figure(figsize=(11, 11))
        sns.heatmap(varmat, square=True, vmin=0, vmax=1)
        plt.show()

    def alignment_plot(self, show_grid=False, save_html=False):
        self._calc_cover_set()

        # Put TR intervals reported by datander on diagonal
        # diagonal
        shapes = [make_line(0, 0, self.read_len, self.read_len, 'grey', 3)]
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
                            textfont=dict(size=10, color="grey"),
                            mode="text",
                            name="interval #")

        data = [trace1, trace2, trace3, trace4]

        if self.encodings is not None and self.read_id in set(self.encodings["read_id"]):
            encodings = self.encodings[self.encodings["read_id"] == self.read_id]

            shapes += list(encodings.apply(lambda df: make_line(df["start"],
                                                                df["start"],
                                                                df["end"],
                                                                df["end"],
                                                                rgb2hex(cm.jet(df["diff"] * 3)),
                                                                5),
                                           axis=1))
            data += [go.Scatter(x=encodings["start"],
                                y=encodings["start"],
                                text=encodings.apply(lambda df: f"{df['peak_id']}:{df['repr_id']}{'+' if df['strand'] == 0 else '-'} ",
                                                     axis=1),
                                textposition="bottom left",
                                textfont=dict(size=10, color="black"),
                                mode="text",
                                name="representative ID"),
                     go.Scatter(x=encodings["start"],
                                y=encodings["start"],
                                text=encodings.apply(lambda df: f"{df['peak_id']}:{df['repr_id']}{'+' if df['strand'] == 0 else '-'}<br>[{df['start']}, {df['end']}] ({df['length']} bp)<br>diff = {df['diff']}",
                                                     axis=1),
                                hoverinfo="text",
                                showlegend=False)]

        layout = go.Layout(width=1000, height=1000,
                           hovermode='closest',
                           xaxis=dict(showgrid=False,
                                      range=[-self.read_len * 0.05, self.read_len + 100]),
                           yaxis=dict(showgrid=False,
                                      range=[0, self.read_len],
                                      autorange='reversed'),
                           shapes=shapes)

        fig = go.Figure(data=data, layout=layout)
        if not save_html:
            py.iplot(fig)
        else:
            py.plot(fig, filename=f"{self.read_id}.alignment_plot.html")

    def path_plot(self, show_snake=True, save_html=False):
        self._load_paths(plot=True, show_snake=show_snake)

        # diagonal
        shapes = [make_line(0, 0, self.read_len, self.read_len, 'yellow', 3)]
        # path trace (and reflecting snakes)
        for path in self.paths:
            shapes.extend(path.shapes)

        # Show alignment paths (and snakes)
        trace1 = go.Scatter(x=[x.ab for x in self.paths],
                            y=[x.bb for x in self.paths],
                            mode='markers',
                            name="start")

        trace2 = go.Scatter(x=[x.ae for x in self.paths],
                            y=[x.be for x in self.paths],
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

        fig = go.Figure(data=[trace1, trace2], layout=layout)
        if not save_html:
            py.iplot(fig)
        else:
            py.plot(fig, filename=f"{self.read_id}.path_plot.html")

    def show_dag(self):
        import networkx as nx
        from networkx.drawing.nx_agraph import graphviz_layout

        self._load_paths()

        for path in self.paths:
            # Shorter than duplication
            if len(path.unit_alignments) == 0:
                continue

            path.unit_consensus()

            print(path.DAG)

            # Show DAG   # TODO: show partial graph
            plt.figure(figsize=(18, 10))
            plt.axis("off")
            #pos = nx.spectral_layout(DAG)
            #pos = nx.circular_layout(DAG)
            #pos = graphviz_layout(DAG, prog="dot")
            pos = graphviz_layout(path.DAG, prog="neato")
            nx.draw_networkx(path.DAG, pos, with_labels=False, node_size=1, font_size=1)   # TODO: output as dot file
            #edge_weights = nx.get_edge_attributes(path.DAG, 'weight')
            #nx.draw_networkx_edge_labels(DAG, pos, edge_labels=edge_weights)
            #plt.show()
