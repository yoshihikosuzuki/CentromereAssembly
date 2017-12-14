import re
import numpy as np
import pandas as pd
from interval import interval
from io import StringIO
from IPython.display import display

from datruf_utils import (run_command,
                          make_line,
                          interval_len,
                          subtract_interval)


# Object datruf is either Viewer or Runner
def load_tr_intervals(datruf):
    command = ("DBdump -mtan %s %d | awk '$1 == \"T0\" {print $0}'"
               % (datruf.db_file, datruf.read_id))
    dbdump = run_command(command).strip().split(' ')
    return [(int(dbdump[1 + 2 * (i + 1)]), int(dbdump[1 + 2 * (i + 1) + 1]))
            for i in range(int(dbdump[1]))]   # [(start_pos, end_pos), ...]


def load_alignments(datruf):
    command = ("LAdump -c %s %s %d | awk '$1 == \"C\" {print $0}'"
               % (datruf.db_file, datruf.las_file, datruf.read_id))
    ladump = StringIO(run_command(command))
    # NOTE: index should be replaced from "C" to interger?
    return pd.read_csv(ladump, sep=" ",
                       names=("abpos", "aepos", "bbpos", "bepos"))

""" TODO: is it better to sort alignments here?
    datruf.alignments = (datruf.alignments
                         .assign(distance=lambda x: x["abpos"] - x["bbpos"])
                         .sort_values(by="abpos", kind="mergesort")
                         .sort_values(by="distance", kind="mergesort"))
"""


def convert_symbol(aseq, bseq, symbols):
    converted = ""
    for i in range(len(symbols)):
        if symbols[i] == '|':   # TODO: replace is faster?
            symbol = "M"
        else:
            if aseq[i] == '-':
                symbol = "I"
            elif bseq[i] == '-':
                symbol = "D"
            else:
                symbol = "N"
        converted += symbol
    return converted


def load_paths(datruf):   # TODO: modify LAshow4pathplot so that only aligned regions are output
    # Load paths of alignments in the minimum cover set
    command = ("LAshow4pathplot -a %s %s %d | sed 's/,//g'"
               "| awk -F'[' 'NF == 1 {print $1} NF == 2 {print $2}'"
               "| awk -F']' '{print $1}'"
               % (datruf.db_file, datruf.las_file, datruf.read_id))
    lashow = run_command(command).strip().split('\n')

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
                #paths.append((ab, ae, bb, be, aseq, bseq, convert_symbol(aseq, bseq, symbols)))
                paths.append(Path(ab, ae, bb, be, Alignment(aseq, bseq, convert_symbol(aseq, bseq, symbols))))
                #print(ab, ae, bb, be, len(aseq))
            ab, ae, bb, be = list(map(int, data[2:6]))
            if (ab, ae, bb, be) in datruf.min_cover_set:
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
        #paths.append((ab, ae, bb, be, aseq, bseq, convert_symbol(aseq, bseq, symbols)))   # TODO: change to dict?
        paths.append(Path(ab, ae, bb, be, Alignment(aseq, bseq, convert_symbol(aseq, bseq, symbols))))

    return paths
