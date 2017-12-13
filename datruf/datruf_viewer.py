import os
import re
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from interval import interval
from io import StringIO

from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib.image as img
import plotly.offline as py
import plotly.graph_objs as go

from datruf_utils import (run_command,
                          make_line,
                          interval_len,
                          subtract_interval)

plt.style.use('ggplot')
py.init_notebook_mode()


class Viewer:
    """
    Dot plot/alignment plot/alignment path plot viewer for one read.
    The 'read_id' is same as the one used in DAZZ_DB.
    This is considered to be used in Jupyter Notebook for now.
    """

    def __init__(self, root_dir, db_file, las_file, out_dir, gepard):
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        self.root_dir = root_dir
        self.db_file = db_file
        self.las_file = las_file
        self.out_dir = out_dir
        self.gepard = gepard

    def set_read(self, read_id):
        self.read_id = read_id

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
            command = ("DBshow %s %s > %s"
                       % (self.db_file, str(self.read_id), out_fasta))
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
        tr_intervals = load_tr_intervals(self.db_file, self.read_id)
        alignments = load_alignments(self.db_file, self.las_file, self.read_id)

        # Put TR intervals reported by datander on diagonal
        # diagonal
        shapes = [make_line(0, 0, self.read_len, self.read_len, 'yellow', 3)]
        for start, end in tr_intervals:
            # TR interval
            shapes.append(make_line(start, start, end, end, 'black', 3))

        alignments, cover_set, additional_shapes = calculate_cover_set(tr_intervals, alignments, self.read_len, show_grid)
        shapes.extend(additional_shapes)
        self.min_cover_set, additional_shapes = calculate_min_cover_set(cover_set)   # TODO: check the existance of min_cover_set first of all
        shapes.extend(additional_shapes)
        
        # Show alignment plot
        trace1 = go.Scatter(x=list(alignments["abpos"]), y=list(alignments["bbpos"]), text=list(range(1, len(alignments) + 1)), textposition="top right", textfont=dict(size=10, color="red"), mode='text', name="alignment #")
        trace2 = go.Scatter(x=list(alignments["abpos"]), y=list(alignments["bbpos"]), mode='markers', name="start")
        trace3 = go.Scatter(x=list(alignments["aepos"]), y=list(alignments["bepos"]), mode='markers', name="end")
        trace4 = go.Scatter(x=[x[0] for x in tr_intervals], y=[x[0] for x in tr_intervals], text=list(range(1, len(tr_intervals) + 1)), textposition="top right", textfont=dict(size=10, color="black"), mode="text", name="interval #")
        layout = go.Layout(width=725, height=725, hovermode='closest', xaxis=dict(showgrid=False, range=[0, self.read_len + 100]), yaxis=dict(autorange='reversed', showgrid=False, range=[0, self.read_len]), shapes=shapes)
        py.iplot(go.Figure(data=[trace1, trace2, trace3, trace4], layout=layout))

    def path_plot(self, read_id):
        # Load paths of alignments in the minimum covering set
        lashow = subprocess.check_output("LAshow4pathplot -a %s %s %s | sed 's/,//g' | awk -F'[' 'NF == 1 {print $1} NF == 2 {print $2}' | awk -F']' '{print $1}'" % (self.db_file, self.las_file, str(read_id)), shell=True).decode('utf-8').strip().split('\n')
        paths = load_paths(lashow, self.min_cover_set)
    
        # Show path plot
        shapes = [{'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': 0, 'y0': 0, 'x1': self.read_len, 'y1': self.read_len, 'line': {'color': 'yellow', 'width': 3}}]   # diagonal
        shapes.extend(generate_path_trace(paths))
        
        trace1 = go.Scatter(x=[x[0] for x in paths], y=[x[2] for x in paths], mode='markers', name="start")
        trace2 = go.Scatter(x=[x[1] for x in paths], y=[x[3] for x in paths], mode='markers', name="end")
        layout = go.Layout(width=1000, height=1000, hovermode='closest', xaxis=dict(showgrid=False, range=[0, self.read_len + 100]), yaxis=dict(autorange='reversed', showgrid=False, range=[0, self.read_len]), shapes=shapes)
        py.iplot(go.Figure(data=[trace1, trace2], layout=layout))
