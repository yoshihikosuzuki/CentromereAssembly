import numpy as np
import pandas as pd
from collections import Counter
from interval import interval
from logzero import logger
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import plotly.offline as py
import plotly.graph_objs as go
from BITS.utils import make_line, interval_len


class Plotter():
    """
    Plot results of tandem repeat analysis. Expected methods are datander, TRF,
    and mTR; however, one can specify arbitrary number of input result files
    using add_result(). The format of input file must be as follows (separater
    is a tab):

        dbid	start	end	unit length	unit sequence
    0	20	5	6897	13	AAGAGAGAAAGAG
    1	53	1	20823	10	AGAATAACAT
    """

    columns = ("dbid",
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
        logger.info(f"end dbid ({name}) = {end_dbid}")
        result = result[result["dbid"] <= end_dbid]
        result = result[result["dbid"] >= self.start_dbid]

        for column in result.columns.values:
            if column not in Plotter.columns:
                result.drop(column, axis=1)

        result["method"] = [name] * len(result)
        result["repeat length"] = result.apply(lambda x: x["end"] - x["start"], axis=1)   # NOTE: be careful of no +1
        result["copy number"] = result.apply(lambda x: x["repeat length"] / x["unit length"], axis=1)

        self.results = (pd.concat([self.results, result])
                        .sort_values(by="start", kind="mergesort")
                        .sort_values(by="dbid", kind="mergesort")
                        .reset_index(drop=True)
                        .loc[:, Plotter.columns])

    def get_results(self):
        return self.results

    def venn_detected_read(self, methods=["datander", "TRF", "mTR"]):
        if len(methods) not in (2, 3):
            return

        reads = [set(self.results[self.results["method"] == method]["dbid"])
                 for method in methods]

        fig, ax = plt.subplots(figsize=(10, 10))
        v = (venn2(reads, methods) if len(methods) == 2
             else venn3(reads, methods))
        for text in v.set_labels:
            text.set_fontsize(15)
        for text in v.subset_labels:
            text.set_fontsize(15)
        ax.set_title("The number of reads in which TR is detected")
        plt.show()

    # Histogram of the frequency of unit length
    def plot_unit_frequency(self):
        methods = sorted(list(set(self.results["method"])))
        data = []

        for i in range(len(methods)):
            result = self.results[self.results["method"] == methods[i]]
            unit_count = Counter()
            for k, v in result.iterrows():
                unit_count[v["unit length"]] += v["copy number"]
            data.append(go.Bar(x=list(unit_count.keys()),
                               y=list(unit_count.values()),
                               opacity=1 if i == 0 else 0.6,
                               name=methods[i]))

        layout = go.Layout(title=methods[0] if len(methods) == 1 else None,
                           xaxis=dict(title="Unit length",
                                      range=[0, 1000]),
                           yaxis=dict(title="Frequency",
                                      range=[0, 2000]))

        py.iplot(go.Figure(data=data, layout=layout))

    # Scatter plot of estimated unit lengths between 2 methods
    def plot_unit_length(self, methods=["datander", "TRF"]):
        if len(methods) != 2:
            return

        def convert_to_list(results, method):
            ret = []
            result = results[results["method"] == method]
            for dbid, intervals in result.groupby("dbid"):
                interval_list = []
                for k, v in intervals.iterrows():
                    interval_list.append(v[["start", "end", "unit length"]])
                ret.append((dbid, interval_list))
            return ret

        m1_list, m2_list = [convert_to_list(self.results, method)
                            for method in methods]

        out = []   # [(unit length by method1, unit length by method2, dbid, (average) interval length)]
        m1_index = 0
        m2_index = 0
        while True:
            if m1_index >= len(m1_list) or m2_index >= len(m2_list):
                break
            if m1_list[m1_index][0] < m2_list[m2_index][0]:
                dbid, m1_data = m1_list[m1_index]
                for m1_interval in m1_data:
                    start, end, unit_len = m1_interval
                    out.append((unit_len, 0, dbid, min(1000, end - start + 1)))
                m1_index += 1
                continue
            if m1_list[m1_index][0] > m2_list[m2_index][0]:
                dbid, m2_data = m2_list[m2_index]
                for m2_interval in m2_data:
                    start, end, unit_len = m2_interval
                    out.append((0, unit_len, dbid, min(1000, end - start + 1)))
                m2_index += 1
                continue

            for m1_interval in m1_list[m1_index][1]:
                for m2_interval in m2_list[m2_index][1]:
                    m1_intvl = interval(m1_interval[:2])
                    m2_intvl = interval(m2_interval[:2])
                if interval_len(m1_intvl & m2_intvl) >= 0.5 * min(interval_len(m1_intvl), interval_len(m2_intvl)):
                    out.append((m1_interval[2],
                                m2_interval[2],
                                m1_list[m1_index][0],
                                min(1000, (interval_len(m1_intvl) + interval_len(m2_intvl)) / 2)))   # TODO: output missing intervals in only one method
            m1_index += 1
            m2_index += 1

        trace = go.Scatter(x=[x[0] for x in out],
                           y=[x[1] for x in out],
                           text=[x[2] for x in out],
                           mode="markers",
                           marker=dict(size=1,
                                       color=[x[3] for x in out],
                                       colorscale='Viridis',
                                       reversescale=True,
                                       showscale=True))

        min_max = min(max([x[0] for x in out]), max([x[1] for x in out]))
        shapes = [make_line(0, 0, min_max, min_max, 'red', 1)]
        layout = go.Layout(title="Estimated unit length [bp] (Color bar = interval length [bp])",
                           xaxis=dict(title=methods[0]),
                           yaxis=dict(title=methods[1]),
                           shapes=shapes,
                           width=950,
                           height=950,
                           hovermode='closest')

        py.iplot(go.Figure(data=[trace], layout=layout))

        # Histogram of the ratio between the estimated unit lengths
        # NOTE: only >50 bp estimated lengths in both methods
        trace = go.Histogram(x=[x[1] / x[0] for x in out
                                if x[0] >= 50 and x[1] >= 50],
                             xbins=dict(start=0, end=5, size=0.01))
        layout = go.Layout(xaxis=dict(title=f"Ratio of length ({methods[1]} / {methods[0]})",
                                      range=[0, 5]),
                           yaxis=dict(title="Frequency"),
                           hovermode='closest',
                           bargap=10,
                           bargroupgap=10)
        py.iplot(go.Figure(data=[trace], layout=layout))
