import igraph as ig
import plotly.graph_objs as go
from BITS.plot.plotly import make_line, show_plot


def overlaps_to_string_graph(overlaps):
    """Construct a string graph from `overlaps` <List[Overlap]>."""
    edges = set()
    for overlap in overlaps:
        f_id, g_id, strand, f_start, f_end, f_len, g_start, g_end, g_len, diff = overlap.astuple()

        # Convert g_[start|end] so that `strand(g.seq[g_start:g_end])` is the overlapping sequence.
        if strand == 1:
            g_start, g_end = g_len - g_end, g_len - g_start

        # Determine the overlap type
        if f_start > 0:
            if strand == 0:
                if g_end == g_len:
                    #    f.B            f.E
                    #  f  -------------->
                    #  g     ------->
                    #      g.B      g.E
                    overlap_type = "contains"
                else:
                    #    f.B       f.E
                    #  f  --------->
                    #  g      ---------->
                    #        g.B        g.E
                    overlap_type = "suffix-prefix"
            else:
                if g_start == 0:
                    #    f.B            f.E
                    #  f  -------------->
                    #  g     <-------
                    #      g.E      g.B
                    overlap_type = "contains"
                else:
                    #    f.B       f.E
                    #  f  --------->
                    #  g      <----------
                    #        g.E        g.B
                    overlap_type = "suffix-suffix"
        else:
            if f_end == f_len:
                overlap_type = "contained"
            else:
                if strand == 0:
                    if g_start == 0:
                        overlap_type = "contains"
                    else:
                        #        f.B        f.E
                        #  f      ---------->
                        #  g  --------->
                        #    g.B       g.E
                        overlap_type = "prefix-suffix"
                else:
                    if g_end == g_len:
                        overlap_type = "contains"
                    else:
                        #        f.B       f.E
                        #  f      --------->
                        #  g  <---------
                        #    g.E       g.B
                        overlap_type = "prefix-prefix"

        # Convert an overlap to nodes and edges in the string graph
        if overlap_type in ["contains", "contained"]:   # contained removal
            continue
        elif overlap_type == "suffix-prefix":
            edges.update([(f"{g_id}:B", f"{f_id}:B", f_start, diff),
                          (f"{f_id}:E", f"{g_id}:E", g_len - g_end, diff)])
        elif overlap_type == "suffix-suffix":
            edges.update([(f"{g_id}:E", f"{f_id}:B", f_start, diff),
                          (f"{f_id}:E", f"{g_id}:B", g_start, diff)])
        elif overlap_type == "prefix-suffix":
            edges.update([(f"{f_id}:B", f"{g_id}:B", g_start, diff),
                          (f"{g_id}:E", f"{f_id}:E", f_len - f_end, diff)])
        else:   # prefix-prefix
            edges.update([(f"{f_id}:B", f"{g_id}:E", g_len - g_end, diff),
                          (f"{g_id}:B", f"{f_id}:E", f_len - f_end, diff)])

    return ig.Graph.DictList(edges=(dict(source=s, target=t, length=l, diff=d) for s, t, l, d in edges),
                             vertices=None,
                             directed=True)


def reduce_transitive_edges(sg, fuzz=100):
    """Perform the transitive edge reduction introduced in [Myers, 2005].
    `fuzz` [bp] determines the degree of fluctuation of overlap start/end positions the algorithm accepts.
    """
    v_mark = ["vacant" for v in sg.vs]
    e_reduce = {e.tuple: False for e in sg.es}

    for v in sg.vs:
        if v.outdegree() == 0:
            continue

        oes = sorted(sg.es.select(_source=v.index), key=lambda x: x["length"])
        longest = oes[-1]["length"] + fuzz
        for oe in oes:
            v_mark[oe.target] = "inplay"

        for oe in oes:
            if v_mark[oe.target] == "inplay":
                ooes = sorted(sg.es.select(_source=oe.target), key=lambda x: x["length"])
                for ooe in ooes:
                    if oe["length"] + ooe["length"] <= longest and v_mark[ooe.target] == "inplay":
                        v_mark[ooe.target] = "eliminated"

        for oe in oes:
            ooes = sorted(sg.es.select(_source=oe.target), key=lambda x: x["length"])
            if len(ooes) > 1:
                shortest = ooes[0].target
                if v_mark[shortest] == "inplay":
                    v_mark[shortest] == "eliminated"
            for ooe in ooes:
                if ooe["length"] < fuzz and v_mark[ooe.target] == "inplay":
                    v_mark[ooe.target] = "eliminated"

        for oe in oes:
            if v_mark[oe.target] == "eliminated":
                e_reduce[oe.tuple] = True
            v_mark[oe.target] = "vacant"

    # Re-construct a graph
    return ig.Graph.DictList(edges=(dict(source=e["source"],
                                         target=e["target"],
                                         length=e["length"],
                                         diff=e["diff"])
                                    for e in sg.es
                                    if not e_reduce[e.tuple]),
                             vertices=None,
                             directed=True)


def draw_string_graph(sg, reads_by_id=None, size=1000):
    # Layout the nodes
    pos = sg.layout('kk')   # {node: (x_coord, y_coord)}

    # Node objects
    cover_rates = None
    if reads_by_id is not None:
        reads = [reads_by_id[int(v["name"].split(':')[0])] for v in sg.vs]
        cover_rates = [sum([unit.length for unit in read.units]) / read.length for read in reads]
        
    trace_node = go.Scatter(x=[pos[v.index][0] for v in sg.vs],
                            y=[pos[v.index][1] for v in sg.vs],
                            text=[f"{v['name']}<br>{round(cover_rates[v.index] * 100, 1) if cover_rates is not None else '-'}% covered"
                                  for v in sg.vs],
                            mode="markers",
                            marker=dict(
                                color=cover_rates,
                                showscale=False,
                                colorscale='YlGnBu',
                                reversescale=False,
                                size=10,
                                line=dict(width=2)))

    # Edge objects
    edges = [e.tuple for e in sg.es]   # List[(source, target)]
    trace_edge = go.Scatter(x=[x for s, t in edges for x in (pos[s][0], pos[t][0], None)],
                            y=[x for s, t in edges for x in (pos[s][1], pos[t][1], None)],
                            mode="lines", line=dict(width=0.5, color='black'))
    shapes = [make_line(pos[s][0] + (pos[t][0] - pos[s][0]) * 0.7,
                        pos[s][1] + (pos[t][1] - pos[s][1]) * 0.7,
                        pos[t][0],
                        pos[t][1],
                        "black", 4, "below")
              for s, t in edges]

    # Draw the graph
    layout = go.Layout(width=size, height=size,
                       xaxis=dict(showgrid=False,
                                  zeroline=False,
                                  showticklabels=True),
                       yaxis=dict(showgrid=False,
                                  zeroline=False,
                                  showticklabels=True),
                       shapes=shapes,
                       hovermode='closest',
                       margin=go.layout.Margin(l=0, r=0, b=0, t=0),
                       showlegend=False)
    show_plot([trace_edge, trace_node], layout)

    return pos
