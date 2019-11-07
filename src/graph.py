import igraph as ig
import plotly.graph_objs as go
from BITS.plot.plotly import make_line, make_scatter, show_plot


def convert_overlap(overlap):
    """Add overlap type after adjusting `g_[start|end]` so that `strand(g.seq[g_start:g_end])`
    is the overlapping sequence.
    """
    f_id, g_id, strand, f_start, f_end, f_len, g_start, g_end, g_len, diff = overlap.astuple()
    assert f_id != g_id, "Self overlap is not allowed"

    if strand == 1:
        g_start, g_end = g_len - g_end, g_len - g_start

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

    return (f_id, g_id, strand, f_start, f_end, f_len, g_start, g_end, g_len, diff, overlap_type)


def overlaps_to_string_graph(overlaps):
    """Construct a string graph from `overlaps` <List[Overlap]>. All the overlaps given here are
    treated as true overlaps. That is, filtering of the overlaps must be finished in advance.
    """
    # Modify `g_[start|end]` and list up contained reads
    converted_overlaps = [convert_overlap(overlap) for overlap in overlaps]
    contained_reads = set()
    for overlap in converted_overlaps:
        if overlap[-1] == "contains":
            contained_reads.add(overlap[1])
        elif overlap[-1] == "contained":
            contained_reads.add(overlap[0])

    # Convert an overlap to nodes and edges in the string graph
    edges = set()
    for overlap in converted_overlaps:
        f_id, g_id, strand, f_start, f_end, f_len, g_start, g_end, g_len, diff, overlap_type = overlap
        
        if f_id in contained_reads or g_id in contained_reads:
            continue

        # TODO: best overlap logic?

        assert overlap_type not in ("contains", "contained"), "Unremoved contained reads"
        if overlap_type == "suffix-prefix":
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


def draw_string_graph(sg, reads_by_id=None, width=1000, height=1000, out_fname=None):
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
    trace_edge = make_scatter(x=[x for e in sg.es for x in (pos[e.source][0], pos[e.target][0], None)],
                              y=[x for e in sg.es for x in (pos[e.source][1], pos[e.target][1], None)],
                              mode="lines", col="black")
    trace_annot = make_scatter(x=[(pos[e.source][0] + pos[e.target][0]) / 2 for e in sg.es],
                               y=[(pos[e.source][1] + pos[e.target][1]) / 2 for e in sg.es],
                               text=[f"{e['length']} bp, {e['diff']}% diff" for e in sg.es],
                               mode="markers", marker_size=1, col="black")
    shapes = [make_line(pos[e.source][0] + (pos[e.target][0] - pos[e.source][0]) * 0.7,
                        pos[e.source][1] + (pos[e.target][1] - pos[e.source][1]) * 0.7,
                        pos[e.target][0],
                        pos[e.target][1],
                        "black",
                        max(e["length"] // 1000, 3),
                        "below")
              for e in sg.es]

    # Draw the graph
    layout = go.Layout(width=width, height=height,
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
    show_plot([trace_edge, trace_annot, trace_node], layout, out_fname=out_fname)

    return pos
