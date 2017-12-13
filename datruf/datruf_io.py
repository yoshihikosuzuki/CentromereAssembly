import os
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


def load_tr_intervals(db_file, read_id):   # TODO: move to another file ("datruf_io.py" or something like that)   # TODO: make 1 read version (load data using read_id) and all read version (already loaded)
    dbdump = subprocess.check_output("DBdump -mtan %s %d | awk '$1 == \"T0\" {print $0}'" % (db_file, read_id), shell=True).decode('utf-8').strip().split(' ')
    return [(int(dbdump[1 + 2 * (i + 1)]), int(dbdump[1 + 2 * (i + 1) + 1])) for i in range(int(dbdump[1]))]   # [(start_pos, end_pos), ...]


def load_alignments(db_file, las_file, read_id):
    ladump = StringIO(subprocess.check_output("LAdump -c %s %s %s | awk '$1 == \"C\" {print $0}'" % (db_file, las_file, str(read_id)), shell=True).decode('utf-8'))#.strip().split('\n'))
    return pd.read_csv(ladump, sep=" ", names=("abpos", "aepos", "bbpos", "bepos"))   # NOTE: index should be replaced from "C" to interger?


def load_paths(lashow, min_cover_set):   # TODO: modify LAshow4pathplot so that only aligned regions are output
    paths = []
    flag_first = True
    for line in lashow:
        data = line.strip().split('\t')
        if len(data) > 1:
            if flag_first:
                flag_first = False
            elif flag_add:
                aseq = aseq[prefix_cut : len(aseq) - suffix_cut]
                bseq = bseq[prefix_cut : len(bseq) - suffix_cut]
                symbols = symbols[:len(symbols) - suffix_cut]
                paths.append((ab, ae, bb, be, aseq, bseq, symbols))
                #print(ab, ae, bb, be, len(aseq))
            ab, ae, bb, be = list(map(int, data[2:6]))
            if (ab, ae, bb, be) in min_cover_set:
                flag_add = True
            else:
                flag_add = False
                continue
            #print(read_id, ab, ae, bb, be)
            aseq = ""
            bseq = ""
            symbols = ""
            counter = 0
            prefix_cut = 0
            suffix_cut = 0
        else:
            if not flag_add:
                continue
            data = data[0]
            if counter % 3 == 0:
                aseq += data
                if '.' in data:
                    if counter == 0:
                        prefix_cut = re.search(r'\.+', data).span()[1]
                    else:
                        suffix_cut += len(data) - re.search(r'\.+', data).span()[0]   # there is a possibility that the "..." region spans 2 rows
            elif counter % 3 == 1:
                symbols += data
                if data[0] == ':':
                    suffix_cut += len(data)   # there is a possibility that the ":::" region spans 2 rows
            else:
                bseq += data
                if '.' in data:
                    if counter == 2:
                        prefix_cut = re.search(r'\.+', data).span()[1]
                    else:
                        suffix_cut += len(data) - re.search(r'\.+', data).span()[0]
            counter += 1
           
    if flag_add:
        aseq = aseq[prefix_cut:len(aseq) - suffix_cut]
        bseq = bseq[prefix_cut:len(bseq) - suffix_cut]
        symbols = symbols[:len(symbols) - suffix_cut]
        paths.append((ab, ae, bb, be, aseq, bseq, symbols))
        
    return paths
