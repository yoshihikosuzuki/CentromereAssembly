import numpy as np
import pandas as pd
import networkx as nx
from interval import interval
from IPython.display import display

from datruf_utils import (make_line,
                          interval_len,
                          subtract_interval)

# -----------------------------------------------------------------------
# TODO: is it better to not use datruf but only data needed as arguments?
# -----------------------------------------------------------------------


def calc_cover_set(datruf):
    datruf.alignments = (datruf.alignments
                         .assign(distance=lambda x: x["abpos"] - x["bbpos"])
                         .sort_values(by="abpos", kind="mergesort")
                         .sort_values(by="distance", kind="mergesort"))   # TODO: should be sorted beforehand?

    uncovered_intervals = interval()
    for intvl in datruf.tr_intervals:
        uncovered_intervals |= interval[intvl]

    cover_set = set()
    for alignment in datruf.alignments.iterrows():
        ab, ae, bb, be = alignment[1][["abpos", "aepos", "bbpos", "bepos"]]
        if not (0.95 <= float(ae - ab) / (be - bb) <= 1.05):   # abnormal slope (TODO: proper?)
            pass
        elif ab - be <= 20:   # TR alignment
            intersect = interval[bb, ae] & uncovered_intervals
            short_intervals = interval()
            for intvl in intersect.components:
                if interval_len(intvl) <= 20:
                    short_intervals |= intvl
            intersect = subtract_interval(intersect, short_intervals, length_threshold=0)   # TODO: meaning?

            # TODO: inspect whether codes below are still proper or not
            flag_deletion = False
            if len(intersect) == 1 and ((bb == intersect[0][0]) != (ae == intersect[-1][1])):
                flag_deletion = True

            if (not flag_deletion and intersect != interval()) or (interval_len(intersect) >= 1.5 * (ab - bb)):   # outer TR, or more than 1 units are in the uncovered regions
                cover_set.add((bb, ae, ab, be))    # NOTE: be careful with the order!
                uncovered_intervals = subtract_interval(uncovered_intervals, interval[bb, ae], length_threshold=100)

    return cover_set


def calc_min_cover_set(cover_set):
    cover_set = sorted(list(cover_set))   # by bbpos
    flanking_length = 100   # for suppressing fluctuation due to the sequencing error
    min_cover_set = set()
    index = 0
    furthest_point = 0
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
            index = next_index + 1
            furthest_point = next_furthest
    return min_cover_set


class Alignment:
    def __init__(self, aseq, bseq, symbol):
        self.aseq = aseq
        self.bseq = bseq
        self.symbol = symbol


class Path:
    """
    Alignments whose paths are to be inspected. They normally belong to
    (minimum) cover set.
    """

    def __init__(self, ab, ae, bb, be, alignment):
        self.ab = ab
        self.ae = ae
        self.bb = bb
        self.be = be
        self.alignment = alignment   # entire alignment

    def split_alignment(self, plot=False, snake=False):
        # split entire alignment into all unit-vs-unit alignments
        ret = trace_alignment(self, plot=plot, snake=snake)
        if plot is True:
            self.unit_alignments, self.unit_len, self.shapes = ret
        else:
            self.unit_alignments, self.unit_len = ret

        # Determine raw unit sequences
        self.unit_seqs = [unit_alignment.bseq.replace('-', '')
                          for unit_alignment in self.unit_alignments]
        self.unit_seqs.append(self.unit_alignments[-1].aseq.replace('-', ''))

    def unit_consensus(self):
        self.DAG = take_consensus(self.unit_alignments)   # TODO: return consensus sequence


# During the trace,
# 1) calculate average/median distance from diagonal
# 2) divide path into unit paths
# 3) generate shapes of the path (if plot is True)
# 4) generate the reflecting snake of the path (if snake is True)
def trace_alignment(path, plot=False, snake=False):
    ab, ae, bb, be, alignment = path.ab, path.ae, path.bb, path.be, path.alignment
    aseq, bseq, symbol = alignment.aseq, alignment.bseq, alignment.symbol

    unit_alignments = []

    # instances of divided data
    unit_aseq = ""
    unit_bseq = ""
    unit_symbol = ""

    apos = ab   # current a-coordinate in the path
    bpos = bb
    distance_list = [ab - bb]   # from diagonal to each point in the path
    reflection_points = [bb, ab]   # equal to border points of units
    reflection_point = ab

    if plot is True:
        shapes = []
        # coordinate of start position of a series of edges of identical type
        # these are used for alignment path plot
        start_apos = ab
        start_bpos = bb

    for i in range(len(symbol)):
        unit_aseq += aseq[i]
        unit_bseq += bseq[i]
        unit_symbol += symbol[i]

        if symbol[i] in ('M', 'N', 'D'):
            apos += 1
        if symbol[i] in ('M', 'N', 'I'):
            bpos += 1

        # If the type of edge will change (or no edge) in the next step,
        # then add a line consisting of continuous edges of same direction
        if plot is True and (i == len(symbol) - 1
                             or abs(ord(symbol[i]) - ord(symbol[i + 1])) > 1):
            shapes.append(make_line(start_apos, start_bpos, apos, bpos, 'black', 1))   # TODO: change color
            start_apos = apos
            start_bpos = bpos

        distance_list.append(apos - bpos)

        # If bpos will step over current reflection point in the next step,
        # then add a reflection point
        if bpos == reflection_point and (i == len(symbol) - 1
                                         or symbol[i + 1] != 'D'):
            unit_alignments.append(Alignment(unit_aseq, unit_bseq, unit_symbol))
            unit_aseq = unit_bseq = unit_symbol = ""
            reflection_points.append(apos)
            reflection_point = apos

    # TODO: remaining partial aseq|bseq|symbol should be added?

    print("ab:", ab, "ae:", ae, "bb:", bb, "be:", be, "ab-bb:", ab - bb,
          "mean:", int(np.mean(distance_list)),
          "median:", int(np.median(distance_list)))

    # Reflecting snake
    if snake is True:
        col = 'green'
        width = 0.2
        rp = reflection_points
        # Initial horizontal line
        shapes.append(make_line(rp[0], rp[0], rp[1], rp[0], col, width))
        for i in range(1, len(rp) - 1):
            # Vertical line
            shapes.append(make_line(rp[i], rp[i - 1], rp[i], rp[i], col, width))
            # Horizontal line
            shapes.append(make_line(rp[i], rp[i], rp[i + 1], rp[i], col, width))
        # Final vertical line
        shapes.append(make_line(rp[-1], rp[-2], rp[-1], rp[-1], col, width))

    ret = [unit_alignments, int(np.mean(distance_list))]
    if plot is True:
        ret.append(shapes)
    return ret


def take_consensus(unit_alignments):

    aseqs = [x.aseq for x in unit_alignments]
    bseqs = [x.bseq for x in unit_alignments]
    symbols = [x.symbol for x in unit_alignments]

    DAG = nx.DiGraph()
    
    # Add backbone
    backbone = bseqs[0].replace('-', '')
    backbone_path = []
    for i in range(len(backbone) - 1):
        DAG.add_edge("%s:%f" % (backbone[i], i + 1), "%s:%f" % (backbone[i + 1], i + 2), weight=1)   # NOTE: *backbone coordinate is 1-index*
        backbone_path.append("%s:%f" % (backbone[i], i + 1))
    backbone_path.append("%s:%f" % (backbone[-1], len(backbone)))
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
                        
            if symb == 'M':
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
    
    return DAG
    # Calculate maximum weighted path
