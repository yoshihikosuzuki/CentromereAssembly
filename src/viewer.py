from dataclasses import dataclass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go
from logzero import logger
from BITS.seq.plot import DotPlot
from BITS.plot.plotly import make_line, make_scatter
from BITS.util.proc import run_command
from .datruf import load_dumps, filter_alignments


@dataclass(repr=False, eq=False)
class Viewer:
    db_file: str
    las_file: str
    out_dir: str
    gepard: str
    encodings: pd.DataFrame = None

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")

    def set_read(self, read_id):
        self.read_id = read_id
        self.read_len = int(run_command(f"DBdump -s {self.db_file} {self.read_id} | awk '$1 == \"S\" {{print $2}}'").strip())

    def show(self,
             read_id,
             dot_plot=False,
             varmat_plot=False,
             alignment_plot=True,
             show_grid=False,   # for alignment plot
             out_html=False):   # instead of drawing in the current cell
        self.set_read(read_id)
        if dot_plot:
            self.dot_plot()
        if varmat_plot:
            self.varmat_plot()
        if alignment_plot:
            self.alignment_plot(show_grid, out_html)

    def alignment_plot(self, show_grid=False, out_html=False):
        # Load TR intervals and self alignments
        assert hasattr(self, "read_id"), "No read ID set"
        tr_intervals, alignments = load_dumps(self.read_id,
                                              self.read_id,
                                              self.db_file,
                                              self.las_file)
        tr_intervals, alignments = tr_intervals[self.read_id], alignments[self.read_id]
        tr_alignments = filter_alignments(tr_intervals, alignments)

        shapes = [make_line(0, 0, self.read_len, self.read_len, 'grey', 3)]   # diagonal (= read)
        for start, end in tr_intervals:   # on diagonal
            shapes.append(make_line(start, start, end, end, 'black', 3))

        for ab, ae, bb, be, distance, slope in alignments:
            if (ab, ae, bb, be) in tr_alignments:
                col, width = 'purple', 3
            else:
                col = 'black' if 0.95 <= slope <= 1.05 else 'yellow'
                width = 1
            shapes.append(make_line(ab, bb, ae, be, col, width))   # self alignment

            if show_grid is True:
                for pos in [(ab, 0, ab, self.read_len),
                            (ae, 0, ae, self.read_len),
                            (0, bb, self.read_len, bb),
                            (0, be, self.read_len, be)]:
                    shapes.append(make_line(*pos, col='green', width=0.2))

        trace1 = go.Scatter(x=[x[0] for x in alignments],   # ab
                            y=[x[2] for x in alignments],   # bb
                            text=list(range(1, len(alignments) + 1)),
                            textposition="top right",
                            textfont=dict(size=10, color="red"),
                            mode='text',
                            name="alignment #")

        trace2 = make_scatter([x[0] for x in alignments],   # ab
                              [x[2] for x in alignments],   # bb
                              name="start")
        
        trace3 = make_scatter([x[1] for x in alignments],   # ae
                              [x[3] for x in alignments],   # be
                              name="end")

        trace4 = go.Scatter(x=[x[0] for x in tr_intervals],
                            y=[x[0] for x in tr_intervals],
                            text=list(range(1, len(tr_intervals) + 1)),
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
                                text=encodings.apply(lambda df: f"{df['peak_id']}:{df['master_id']}:{df['repr_id']}{'+' if df['strand'] == 0 else '-'} ",
                                                     axis=1),
                                textposition="bottom left",
                                textfont=dict(size=10, color="black"),
                                mode="text",
                                name="representative ID"),
                     go.Scatter(x=encodings["start"],
                                y=encodings["start"],
                                text=encodings.reset_index().apply(lambda df: f"{df.name} ({df['peak_id']}:{df['master_id']}:{df['repr_id']}{'+' if df['strand'] == 0 else '-'})<br>[{df['start']}, {df['end']}] ({df['length']} bp)<br>diff = {df['diff']}",
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
        if not out_html:
            py.iplot(fig)
        else:
            py.plot(fig, filename=f"{self.read_id}.alignment_plot.html")

    def dot_plot(self):
        out_fasta = f"{self.out_dir}/{self.read_id}.fasta"
        run_command(f"DBshow {self.db_file} {self.read_id} > {out_fasta}")
        DotPlot(self.out_dir, self.gepard).plot_fasta(out_fasta, out_fasta)

    def varmat_plot(self, colname="var_vec_global0.0"):   # TODO: change to integrated var. vec. viewer
        assert self.encodings is None, "Encodings must be specified"
        read_df = self.encodings.pipe(lambda df: df[df["read_id"] == self.read_id]) \
                                .pipe(lambda df: df[df["type"] == "complete"])
        if read_df.shape[0] == 0:
            logger.info("The read does not have any complete variant vector unit. Exit.")
            return
        varmat = np.array(list(read_df[colname].values))
        plt.figure(figsize=(11, 11))
        sns.heatmap(varmat, square=True, vmin=0, vmax=1)
        plt.show()
