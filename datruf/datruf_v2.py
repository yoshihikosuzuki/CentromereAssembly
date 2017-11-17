#!/usr/bin/env python
#coding: utf-8

#import os
import re
import math
#import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict
from interval import interval

#from IPython.display import display
#import matplotlib.pyplot as plt
#import matplotlib.image as img
#from matplotlib_venn import venn2, venn3
#plt.style.use('ggplot')
#import plotly.offline as py
#py.init_notebook_mode()
#import plotly.figure_factory as ff
#import plotly.graph_objs as go


## all datasets
datander_dbdump = "datander_dbdump"   # datander intervals for all reads
datander_ladump = "datander_ladump"   # datander alignments for all reads
#trf_all_outfile = "reads_all.1.1.2.80.10.200.2000.dat"   # TRF result against all reads
#mtr_outfile = "mtr_out"   # mTR results against all reads

## outputs (TR intervals and unit length) by each method
datander_result_fname = "datander_result"
#trf_result_fname = "trf_result"
#mtr_result_fname = "mtr_result"


## convertors between DB ids and fasta headers
dbid_header = pd.read_csv("dbid_header", sep="\t", header=None)
dbid_to_header = dict(zip(map(str, list(dbid_header.iloc[:, 0])), list(dbid_header.iloc[:,1])))
header_to_dbid = {v:k for k, v in dbid_to_header.items()}


#max_dbid = 10000   # for using subset (Be careful it is *int* type)
max_dbid = max(dbid_header.loc[:, 0])   # for using all the datasets

## universal utilities

def interval_len(intvl):
    return intvl[0][1] - intvl[0][0] + 1 if intvl != interval() else 0

# subtraction of integer intervals: A - B, and afterwards remove intervals whose length is less than <length threshold>
def subtract_interval(a_interval, b_interval, length_threshold=0):
    ret_interval = interval()
    intersection = a_interval & b_interval
    a_index = 0
    i_index = 0
    flag_load_new_interval = True   # use when there might be other intersecting intervals in an A-interval
    while a_index < len(a_interval) and i_index < len(intersection):
        if flag_load_new_interval:
            a_intvl = a_interval[a_index]
        else:
            flag_load_new_interval = False
        #print(a_index, i_index)
        #print(a_intvl, intersection[i_index])

        if a_intvl[1] < intersection[i_index][0]:
            ret_interval |= interval(a_intvl)
            a_index += 1
        elif a_intvl[0] > intersection[i_index][1]:
            i_index += 1
        else:
            a_start, a_end = a_intvl
            i_start, i_end = intersection[i_index]
            start_contained = True if min(a_start, i_start) == i_start else False
            end_contained = True if max(a_end, i_end) == i_end else False
            #print(start_contained, end_contained)
            if start_contained:
                if end_contained:
                    a_index += 1
                    i_index += 1
                else:
                    a_intvl = interval[i_end + 1, a_end]   # the tail interval that did not intersect with this B-interval (but maybe with the next B-interval)
                    flag_load_new_interval = False
                    i_index += 1
            else:
                if end_contained:
                    ret_interval |= interval[a_start, i_start - 1]
                    a_index += 1
                else:
                    ret_interval |= interval[a_start, i_start - 1]
                    a_intvl = interval[i_end + 1, a_end]
                    flag_load_new_interval = False
                    i_index += 1
    if not flag_load_new_interval:   # add a tail interval if it exists
        ret_interval |= interval(a_intvl)
        a_index += 1
    while a_index < len(a_interval):   # add remaining intervals in A
        ret_interval |= interval(a_interval[a_index])
        a_index += 1
    ret_interval_long = interval()   # length threshold
    for intvl in ret_interval.components:
        if interval_len(intvl) >= length_threshold:
            ret_interval_long |= intvl
    return ret_interval_long

# for adding an element (i.e. a row) into a dictionary that will be in the end converted to pd.DataFrame
def add_element(dic, columns, counter, data):
    for i in range(len(columns)):
        dic[columns[i]][counter] = data[i]


## load self-vs-self alignments by datander

datander_alignments = defaultdict(dict)
datander_alignments_columns = ("dbid", "abpos", "aepos", "bbpos", "bepos")
counter = 0
with open(datander_ladump, 'r') as f:
    for line in f:
        data = line.strip().split(' ')
        if data[0] == "P":
            dbid = data[1]
            if int(dbid) > max_dbid:   # for using subset
                break
        elif data[0] == "C":
            add_element(datander_alignments, datander_alignments_columns, counter, [dbid] + list(map(int, data[1:5])))
            counter += 1
            
datander_alignments = pd.DataFrame.from_dict(datander_alignments)
datander_alignments = datander_alignments.loc[:, datander_alignments_columns]


## load TR intervals determined by TANmask

datander_intervals = defaultdict(dict)
datander_intervals_columns = ("dbid", "start", "end")
counter = 0
with open(datander_dbdump, 'r') as f:
    for line in f:
        data = line.strip().split(' ')
        if data[0] == "R":
            dbid = data[1]
            if int(dbid) > max_dbid:   # for using subset
                break
        elif data[0] == "T0" and int(data[1]) > 0:
            for i in range(int(data[1])):
                add_element(datander_intervals, datander_intervals_columns, counter, [dbid, int(data[1 + 2 * (i + 1)]), int(data[1 + 2 * (i + 1) + 1])])
                counter += 1
                
datander_intervals = pd.DataFrame.from_dict(datander_intervals)
datander_intervals = datander_intervals.loc[:, datander_intervals_columns]


datander_results = defaultdict(dict)
datander_results_columns = ("dbid", "header", "start", "end", "unit length")
df_counter = 0

for dbid, alignments in datander_alignments.groupby("dbid"):
    if len(datander_intervals[datander_intervals["dbid"] == dbid]) == 0:   # no TR intervals
        continue
    uncovered_intervals = interval()
    for intvl in list(datander_intervals[datander_intervals["dbid"] == dbid].apply(lambda x: (x["start"], x["end"]), axis = 1)):
        uncovered_intervals |= interval[intvl]
    
    alignments = alignments.assign(distance = lambda x: x["abpos"] - x["bbpos"]).sort_values(by = "abpos", kind = "mergesort").sort_values(by = "distance", kind = "mergesort")
    #display(alignments)
    
    #print(header)
    #print(uncovered_intervals)
    
    # collect TR alignments that can cover the TR intervals
    covered_interval = interval()
    cover_set = []
    for alignment in alignments.iterrows():
        ab, ae, bb, be = alignment[1][1:5]
        if ab - be <= 20:
            if interval[bb, ae] & uncovered_intervals != interval():
                cover_set.append((bb, ae, ab, be))    # NOTE: be careful with the order!
                uncovered_intervals = subtract_interval(uncovered_intervals, interval[bb, ae], length_threshold=500)
                if uncovered_intervals == interval():
                    break
                #print("add:", interval[bb, ae])
                #print(uncovered_intervals)
    cover_set = sorted(cover_set)
    #print(cover_set)
    
    # determine minimum alignment cover set
    # contained TR alignmetns are also removed here (TODO: is it really desirable?)
    flanking_length = 100   # for suppressing fluctuation due to the sequencing error
    min_cover_set = []
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
            min_cover_set.append(cover_set[next_index])
            index = next_index + 1
            furthest_point = next_furthest
    #print(min_cover_set)
    
    # calculate unit length from each TR alignment itself in the minimum cover set
    for pos in min_cover_set:
        bb, ae, ab, be = pos
        add_element(datander_results, datander_results_columns, df_counter, (int(dbid), dbid_to_header[dbid], bb, ae, ab - bb))
        df_counter += 1
    
datander_results = pd.DataFrame.from_dict(datander_results)
datander_results = datander_results.loc[:, datander_results_columns]
datander_results = datander_results.sort_values(by="dbid").reset_index(drop=True)


datander_results.to_csv(datander_result_fname, sep="\t")
