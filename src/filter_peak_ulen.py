from logzero import logger
from interval import interval
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
import plotly.offline as py
import plotly.graph_objs as go
from BITS.utils import run_command


def load_db(db_file, tr_read_ids):
    # NOTE: header must not contain tab
    return pd.DataFrame.from_dict({dbid: (*line.split('\t'), tr_read_ids[dbid])   # header, sequence, peak_id
                                   for dbid, line in enumerate(run_command(f"DBshow {db_file} | awk 'BEGIN {{first = 1}} {{if (substr($0, 1, 1) == \">\") {{if (first == 1) {{first = 0}} else {{printf(\"%s\\t%s\\n\", header, seq)}}; header = substr($0, 2); seq = \"\";}} else {{seq = seq $0}}}} END {{printf(\"%s\\t%s\\n\", header, seq)}}'").strip().split('\n'), start=1)
                                   if dbid in tr_read_ids.keys()},
                                  orient="index",
                                  columns=("header", "sequence", "peak_id"))


def find_peak_ulen(tr_file,   # output of datruf
                   db_file,
                   min_ulen=50,
                   max_ulen=500,
                   min_cover_rate=0.8,   # of read by all TRs
                   band_width=5,   # for KDE
                   min_density=0.005,   # of peaks in KDE
                   deviation=0.08,   # <peak_ulen> * (1 +- <deviation>) will be each peak interval
                   plot=False):
    # Detect peak unit lengths with KDE
    read_lengths = {dbid: int(read_length) for dbid, read_length in enumerate(run_command(f"DBdump -s {db_file} | awk '$1 == \"S\" {{print $2}}'").strip().split('\n'), start=1)}
    trs = pd.read_csv(tr_file, sep='\t', index_col=0)
    unit_lens = np.array([end - start
                          for read_id, df in trs.groupby("read_id")
                          if (df["end"] - df["start"]).sum() >= min_cover_rate * read_lengths[read_id]
                          for units in df["units"]
                          for start, end in eval(units)
                          if min_ulen <= end - start <= max_ulen])
    dens = np.exp(KernelDensity(kernel='gaussian', bandwidth=band_width)
                  .fit(unit_lens.reshape(-1, 1))
                  .score_samples(np.arange(min_ulen, max_ulen + 1).reshape(-1, 1)))
    if plot:   # TODO: show peak intervals below
        original = go.Histogram(x=unit_lens,
                                xbins=dict(start=min_ulen, end=max_ulen, size=1))
        layout = go.Layout(xaxis=dict(title='Unit length'),
                           yaxis=dict(title='Frequency'))
        py.iplot(go.Figure(data=[original], layout=layout))
        smoothed = go.Scatter(x=np.arange(min_ulen, max_ulen + 1),
                              y=dens)
        layout = go.Layout(xaxis=dict(title='Unit length'),
                           yaxis=dict(title='Density'))
        py.iplot(go.Figure(data=[smoothed], layout=layout))
        return

    # Determine intervals of unit length in each of which TRs will be analyzed further
    peak_ulens = [min_ulen + i for i in range(1, max_ulen - min_ulen)
                  if dens[i] > dens[i - 1] and dens[i] > dens[i + 1] and dens[i] >= min_density]
    logger.info("Peak unit lengths: " + ', '.join([f"{peak_ulen} bp" for peak_ulen in peak_ulens]))
    peak_intvls = interval(*[[-(- peak_ulen * (1. - deviation) // 1),
                              int(peak_ulen * (1. + deviation))]
                             for peak_ulen in peak_ulens])
    logger.info("Peak intervals: " + ', '.join([f"{start}-{end} bp"
                                                for peak_intvl in peak_intvls.components
                                                for start, end in peak_intvl]))

    # Filter reads for each peak interval
    # NOTE: every read can belong to at most 1 peak class if <min_coevr_rate> is bigger than 0.5
    tr_read_ids = {peak_id: {read_id for read_id, df in trs[trs.apply(lambda df: df["mean_ulen"] in peak_intvl, axis=1)].groupby("read_id")
                             if (df["end"] - df["start"]).sum() >= min_cover_rate * read_lengths[read_id]}
                   for peak_id, peak_intvl in enumerate(peak_intvls.components, start=1)}
    tr_read_ids = {read_id: peak_id
                   for peak_id, read_ids in tr_read_ids.items()
                   for read_id in read_ids}
    tr_reads = load_db(db_file, tr_read_ids)
    tr_reads.to_csv("tr_reads", sep='\t')
