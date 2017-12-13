import re
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from interval import interval
from io import StringIO
from IPython.display import display

from datruf_utils import (run_command,
                          make_line,
                          interval_len,
                          subtract_interval)


def calculate_cover_set(tr_intervals, alignments, read_len, show_grid):
    alignments = alignments.assign(distance=lambda x: x["abpos"] - x["bbpos"]).sort_values(by="abpos", kind="mergesort").sort_values(by="distance", kind="mergesort")
    uncovered_intervals = interval()
    for intvl in tr_intervals:
        uncovered_intervals |= interval[intvl]
    cover_set = []
    shapes = []
    for alignment in alignments.iterrows():
        ab, ae, bb, be = alignment[1][["abpos", "aepos", "bbpos", "bepos"]]
        if not (0.95 <= float(ae - ab) / (be - bb) <= 1.05):   # abnormal slope
            color = 'yellow'
        elif ab - be <= 20:   # TR alignment
            intersect = interval[bb, ae] & uncovered_intervals
            short_intervals = interval()
            for intvl in intersect.components:
                if interval_len(intvl) <= 20:
                    short_intervals |= intvl
            intersect = subtract_interval(intersect, short_intervals, length_threshold=0)

            # TODO: inspect whether codes below are OK or not
            flag_deletion = False
            if len(intersect) == 1 and ((bb == intersect[0][0]) != (ae == intersect[-1][1])):
                flag_deletion = True
                
            #display(alignment)
            #print(uncovered_intervals)
            #print(intersect)
            #print(flag_deletion)
            
            if (not flag_deletion and intersect != interval()) or (interval_len(intersect) >= 1.5 * (ab - bb)):   # outer TR, or more than 1 units are in the uncovered regions
                cover_set.append((bb, ae, ab, be))    # NOTE: be careful with the order!
                uncovered_intervals = subtract_interval(uncovered_intervals, interval[bb, ae], length_threshold=100)
                color = 'red'   # TR alignment and cover set
            else:
                color = 'blue'   # TR alignment but not cover set
        else:
            color = 'black'   # not TR alignment
            
        #print(color)
        shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': ab, 'y0': bb, 'x1': ae, 'y1': be, 'line': {'width': 1, 'color': color}})   # alignments
        
        if show_grid is True:
            line_width = 0.2
            shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': ab, 'y0': 0, 'x1': ab, 'y1': read_len, 'line': {'width': line_width, 'color': 'green'}})
            shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': ae, 'y0': 0, 'x1': ae, 'y1': read_len, 'line': {'width': line_width, 'color': 'green'}})
            shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': 0, 'y0': bb, 'x1': read_len, 'y1': bb, 'line': {'width': line_width, 'color': 'green'}})
            shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': 0, 'y0': be, 'x1': read_len, 'y1': be, 'line': {'width': line_width, 'color': 'green'}})

    return (alignments, cover_set, shapes)


def calculate_min_cover_set(cover_set):
    cover_set = sorted(cover_set)   # by bbpos
    flanking_length = 100   # for suppressing fluctuation due to the sequencing error
    min_cover_set = set()
    index = 0
    furthest_point = 0
    shapes = []
    while index < len(cover_set):
        next_index = -1
        next_furthest = furthest_point
        skip_index = -1
        for i in range(index, len(cover_set)):
            bb, ae, ab, be = cover_set[i]
            if furthest_point < bb:
                skip_index = i
                skip_furthest = bb
                break
            else:
                if next_furthest + flanking_length < ae:
                    next_index = i
                    next_furthest = ae
        if next_index == -1:
            if skip_index == -1:
                break
            index = skip_index
            furthest_point = skip_furthest
        else:
            bb, ae, ab, be = cover_set[next_index]
            min_cover_set.add((ab, ae, bb, be))
            shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': ab, 'y0': bb, 'x1': ae, 'y1': be, 'line': {'width': 3, 'color': 'purple'}})   # minimum cover set
            index = next_index + 1
            furthest_point = next_furthest
    min_cover_set = sorted(list(min_cover_set))
    #print(min_cover_set)
    return (min_cover_set, shapes)


def take_consensus(aseqs, bseqs, symbols):
    DAG = nx.DiGraph()
    
    # Add backbone
    backbone = bseqs[0].replace('-', '')
    backbone_path = []
    for i in range(len(backbone) - 1):
        DAG.add_edge("%s:%f" % (backbone[i], i + 1), "%s:%f" % (backbone[i + 1], i + 2), weight=1)   # NOTE: *backbone coordinate is 1-index*
        backbone_path.append("%s:%f" % (backbone[i], i + 1))
    backbone_path.append("%s:%f" % (backbone[i + 1], i + 2))
    #print("backbone_path:", backbone_path)
    last_coordinate = len(backbone) + 1   # [0, last_coordinate]　is used for node coordinates ([1, last_coordinate - 1] is backbone's interval)
    
    # Iteratively add alignments
    for i in range(len(symbols)):
        bseq = bseqs[i]   # backbone (+ gaps)
        aseq = aseqs[i]   # to be added
        symbol = symbols[i]
        
        #print("bseq_len:", len(bseq), "b_bases", len(bseq.replace('-', '')))
        #print(bseq[-20:])
        #print("aseq_len:", len(aseq), "a_bases", len(aseq.replace('-', '')))
        #print(aseq[-20:])
        #print("symbol_len:", len(symbol))
        #print(symbol[-20:])
        
        b_index = 0
        new_path = []   # will be the next backbone
        bases_between_matches = ""   # coordinates of nodes between two nearest matches between A and B are equally divided by # of the bases
        branch_start_node = ""   # must be a node of backbone. this is also used as flag of first node
        for k in range(len(symbol)):
            symb = symbol[k]
            
            #print(bseq[k], symbol[k], aseq[k])
            #print("b_index:", b_index)
                        
            if symb == '|':
                if branch_start_node == "":   # first match
                    if len(bases_between_matches) == 0:   # first base is first match
                        pass
                        #branch_start_node = backbone_path[b_index]
                    else:   # add first series of I/MM
                        branch_start_pos = 0
                        branch_end_pos = float(backbone_path[b_index].split(':')[1])
                        step_size = (branch_end_pos - branch_start_pos) / (len(bases_between_matches) + 1)
                        source_node = "%s:%f" % (bases_between_matches[0], branch_start_pos + step_size)
                        new_path.append(source_node)
                        for j in range(1, len(bases_between_matches)):
                            target_node = "%s:%f" % (bases_between_matches[j], branch_start_pos + step_size * (j + 1))
                            if DAG.has_edge(source_node, target_node):   # TODO: make function
                                DAG[source_node][target_node]['weight'] += 1
                            else:
                                DAG.add_edge(source_node, target_node, weight=1)
                            new_path.append(target_node)
                            source_node = target_node
                        if DAG.has_edge(source_node, backbone_path[b_index]):
                            DAG[source_node][backbone_path[b_index]]['weight'] += 1
                        else:
                            DAG.add_edge(source_node, backbone_path[b_index], weight=1)
                        #new_path.append(backbone_path[b_index])
                else:   # previous match exists
                    if len(bases_between_matches) == 0:   # no I/MM (that is, continuous matches or deletions)
                        if DAG.has_edge(branch_start_node, backbone_path[b_index]):
                            DAG[branch_start_node][backbone_path[b_index]]['weight'] += 1
                        else:
                            DAG.add_edge(branch_start_node, backbone_path[b_index], weight=1)
                        #new_path.append(backbone_path[b_index])
                    else:   # add a series of I/MM
                        branch_start_pos = float(branch_start_node.split(':')[1])
                        branch_end_pos = float(backbone_path[b_index].split(':')[1])
                        step_size = (branch_end_pos - branch_start_pos) / (len(bases_between_matches) + 1)
                        source_node = branch_start_node
                        for j in range(len(bases_between_matches)):
                            target_node = "%s:%f" % (bases_between_matches[j], branch_start_pos + step_size * (j + 1))
                            if DAG.has_edge(source_node, target_node):
                                DAG[source_node][target_node]['weight'] += 1
                            else:
                                DAG.add_edge(source_node, target_node, weight=1)
                            new_path.append(target_node)
                            source_node = target_node
                        if DAG.has_edge(source_node, backbone_path[b_index]):
                            DAG[source_node][backbone_path[b_index]]['weight'] += 1
                        else:
                            DAG.add_edge(source_node, backbone_path[b_index], weight=1)
                new_path.append(backbone_path[b_index])
                branch_start_node = backbone_path[b_index]
                bases_between_matches = ""
                b_index += 1
            else:
                abase = aseq[k]
                bbase = bseq[k]
                if abase != '-':   # insertion or mismatch
                    bases_between_matches += aseq[k]
                if bbase != '-':
                    b_index += 1
                    
            #print("branch_start_node:", branch_start_node)
            #print("bases_between_matches:", bases_between_matches)
            #print("new_path:", new_path)
            #print("---")
            
        #print("remaining_bases:", bases_between_matches)
        
        # add the last series of I/MM
        if branch_start_node == "":   # no matches between A and B (rare; unit length must be very short)
            branch_start_pos = 0
            branch_end_pos = last_coordinate
            step_size = (branch_end_pos - branch_start_pos) / (len(bases_between_matches) + 1)
            source_node = "%s:%f" % (bases_between_matches[0], branch_start_pos + step_size)
            new_path.append(source_node)
            for j in range(1, len(bases_between_matches)):
                target_node = "%s:%f" % (bases_between_matches[j], branch_start_pos + step_size * (j + 1))
                if DAG.has_edge(source_node, target_node):
                    DAG[source_node][target_node]['weight'] += 1
                else:
                    DAG.add_edge(source_node, target_node, weight=1)
                new_path.append(target_node)
                source_node = target_node
        elif len(bases_between_matches) != 0:   # I/MM exist after previous match
            branch_start_pos = float(branch_start_node.split(':')[1])
            branch_end_pos = last_coordinate
            step_size = (branch_end_pos - branch_start_pos) / (len(bases_between_matches) + 1)
            #print("last_branch_start:", branch_start_pos)
            #print("last_branch_end:", branch_end_pos)
            #print("step_size:", step_size)
            source_node = branch_start_node
            for j in range(len(bases_between_matches)):
                target_node = "%s:%f" % (bases_between_matches[j], branch_start_pos + step_size * (j + 1))
                if DAG.has_edge(source_node, target_node):
                    DAG[source_node][target_node]['weight'] += 1
                else:
                    DAG.add_edge(source_node, target_node, weight=1)
                new_path.append(target_node)
                source_node = target_node
                
        #print("new_path_len (= a_bases):", len(new_path))
        #print(new_path[-20:])
        #print("---")
        backbone_path = new_path
    
    # Draw
    plt.figure(figsize=(18,10))
    plt.axis("off")
    #pos = nx.spectral_layout(DAG)
    #pos = nx.circular_layout(DAG)
    #pos = graphviz_layout(DAG, prog="dot")
    pos = graphviz_layout(DAG, prog="neato")
    edge_weights = nx.get_edge_attributes(DAG, 'weight')
    nx.draw_networkx(DAG, pos, with_labels=False, node_size=1, font_size=1)   # TODO: output as dot file
    #nx.draw_networkx_edge_labels(DAG, pos, edge_labels=edge_weights)
    
    # Calculate maximum weighted path


def generate_path_trace(paths):   # and also calculate average distance(s) (and also take consensu of the units) TODO: divide the functions
    shapes = []
    for path in paths:
        aseqs = []
        bseqs = []
        symbols = []
        unit_aseq = ""
        unit_bseq = ""
        unit_symbol = ""
        
        ab, ae, bb, be, aseq, bseq, symbol = path
        aend = ab
        bend = bb
        astart = aend
        bstart = bend
        distance_list = [aend - bend]
        unit_borders = [bb, ab]
        reflection_point = ab
        prev_edge = -1   # 1 for match/mismatch, 0 for in/del
        for i in range(len(aseq)):
            unit_aseq += aseq[i]
            unit_bseq += bseq[i]
            unit_symbol += symbol[i]
            
            if symbol[i] == '|':
                aend += 1
                bend += 1
                distance_list.append(aend - bend)
                flag_prev_match = 1
                prev_edge = 1
                if i == len(aseq) - 1:
                    shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': astart, 'y0': bstart, 'x1': aend, 'y1': bend, 'line': {'color': 'black', 'width': 1}})
            else:
                if flag_prev_match == 1:
                    shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': astart, 'y0': bstart, 'x1': aend, 'y1': bend, 'line': {'color': 'black', 'width': 1}})
                    flag_prev_match = 0
                astart = aend
                bstart = bend
                if aseq[i] == '-':
                    bend += 1
                    prev_edge = 0
                    col = 'blue'
                elif bseq[i] == '-':
                    aend += 1
                    prev_edge = 0
                    col = 'blue'
                else:
                    aend += 1
                    bend += 1
                    prev_edge = 1
                    col = 'red'
                shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': astart, 'y0': bstart, 'x1': aend, 'y1': bend, 'line': {'color': col, 'width': 1}})
                distance_list.append(aend - bend)
                astart = aend
                bstart = bend
            if bend > reflection_point:
                aseqs.append(unit_aseq[:-1])
                bseqs.append(unit_bseq[:-1])
                symbols.append(unit_symbol[:-1])
                unit_aseq = aseq[i]
                unit_bseq = bseq[i]
                unit_symbol = symbol[i]
                
                if prev_edge == 1:
                    unit_borders.append(aend - 1)
                    reflection_point = aend - 1
                elif prev_edge == 0:
                    unit_borders.append(aend)
                    reflection_point = aend
        print("ab:", ab, "ae:", ae, "bb:", bb, "be:", be, "ab-bb:", ab - bb, "mean:", int(np.mean(distance_list)), "median:", int(np.median(distance_list)))

        # Reflecting snake
        #print(unit_borders)

        shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': unit_borders[0], 'y0': unit_borders[0], 'x1': unit_borders[1], 'y1': unit_borders[0], 'line': {'width': 0.2, 'color': 'green'}})
        for i in range(1, len(unit_borders) - 1):
            shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': unit_borders[i], 'y0': unit_borders[i - 1], 'x1': unit_borders[i], 'y1': unit_borders[i], 'line': {'width': 0.2, 'color': 'green'}})
            shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': unit_borders[i], 'y0': unit_borders[i], 'x1': unit_borders[i + 1], 'y1': unit_borders[i], 'line': {'width': 0.2, 'color': 'green'}})
        shapes.append({'type': 'line', 'xref': 'x', 'yref': 'y', 'x0': unit_borders[-1], 'y0': unit_borders[-2], 'x1': unit_borders[-1], 'y1': unit_borders[-1], 'line': {'width': 0.2, 'color': 'green'}})
        
        # Cut out modules (units) and take consensus
        #take_consensus(aseqs, bseqs, symbols)#, unit_borders)

    return shapes
