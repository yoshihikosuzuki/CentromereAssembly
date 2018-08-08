import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go

from .damaster_core import Runner

plt.style.use('ggplot')


class PeaksViewer:
    def __init__(self, unit_fasta):
        self.runner = Runner(unit_fasta, None, None)
        self.runner.smooth_unit_len_dist()   # NOTE: do self.runner.run() if you want peaks as well

    def plot_unit_len_dist(self):
        # Raw unit length distribution
        data = [go.Histogram(x=self.runner.unit_lens,
                             xbins=dict(start=self.runner.min_len,
                                        end=self.runner.max_len,
                                        size=1))]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Frequency"))
        py.iplot(go.Figure(data=data, layout=layout))

        # Smoothed distribution
        data = [go.Scatter(x=self.runner.ls, y=self.runner.dens)]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Density"))
        py.iplot(go.Figure(data=data, layout=layout))
