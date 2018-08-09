import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go

from .damaster_core import PeaksFinder

plt.style.use('ggplot')


class PeaksViewer:
    def __init__(self, unit_fasta):
        self.finder = PeaksFinder(unit_fasta)
        self.finder.smooth_unit_len_dist()
        # NOTE: do self.finder.run() if you want peak lengths as well

    def plot_unit_len_dist(self):
        # Raw unit length distribution
        data = [go.Histogram(x=self.finder.unit_lens,
                             xbins=dict(start=self.finder.min_len,
                                        end=self.finder.max_len,
                                        size=1))]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Frequency"))
        py.iplot(go.Figure(data=data, layout=layout))

        # Smoothed distribution
        data = [go.Scatter(x=self.finder.ls, y=self.finder.dens)]
        layout = go.Layout(xaxis=dict(title="Unit length"),
                           yaxis=dict(title="Density"))
        py.iplot(go.Figure(data=data, layout=layout))
