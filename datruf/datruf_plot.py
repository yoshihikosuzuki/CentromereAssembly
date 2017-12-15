import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib.image as img
from matplotlib_venn import venn2, venn3
import plotly.offline as py
import plotly.graph_objs as go

plt.style.use('ggplot')
py.init_notebook_mode()


class Plotter():
    """
    Plot results of tandem repeat analysis. Expected methods are datander, TRF,
    and mTR; however, one can specify arbitrary number of input result files
    using add_result(). The format of input file must be as follows (separater
    is a tab):

        dbid	header	start	end	unit length	unit sequence
    0	20	L416/2699/0_6897 RQ=0.826	5	6897	13	AAGAGAGAAAGAG
    1	53	L416/4303/0_23522 RQ=0.866	1	20823	10	AGAATAACAT
    """

    columns = ("dbid",
               "header",   # TODO: remove header?
               "method",
               "start",
               "end",
               "repeat length",
               "unit length",
               "copy number",
               "unit sequence")

    def __init__(self, start_dbid=1, end_dbid=-1):
        self.start_dbid = start_dbid
        self.end_dbid = end_dbid
        self.results = pd.DataFrame(columns=Plotter.columns)

    def add_result(self, name, fname):
        result = pd.read_csv(fname, sep="\t", index_col=0)
        end_dbid = max(result["dbid"]) if self.end_dbid < 1 else self.end_dbid
        print("end dbid =", end_dbid)
        result = result[result["dbid"] <= end_dbid]
        result = result[result["dbid"] >= self.start_dbid]

        result["method"] = [name] * len(result)
        result["repeat length"] = result.apply(lambda x: x["end"] - x["start"], axis=1)   # NOTE: be careful of no +1
        result["copy number"] = result.apply(lambda x: x["repeat length"] / x["unit length"], axis=1)

        self.results = (pd.concat([self.results, result])
                        .sort_values(by="dbid")
                        .reset_index(drop=True)
                        .loc[:, Plotter.columns])

    def show_results(self):
        display(self.results)

    # Histogram of the frequency of unit length
    def plot_unit_frequency(self):
        names = sorted(list(set(self.results["method"])))
        data = []

        for i in range(len(names)):
            result = self.results[self.results["method"] == names[i]]
            unit_count = Counter()
            for k, v in result.iterrows():
                unit_count[v["unit length"]] += v["copy number"]
            data.append(go.Bar(x=list(unit_count.keys()),
                               y=list(unit_count.values()),
                               opacity=1 if i == 0 else 0.6,
                               name=names[i]))

        layout = go.Layout(title=names[0] if len(names) == 1 else None,
                           xaxis=dict(title="Unit length",
                                      range=[0, 1000]),
                           yaxis=dict(title="Frequency",
                                      range=[0, 2000]))

        py.iplot(go.Figure(data=data, layout=layout))
