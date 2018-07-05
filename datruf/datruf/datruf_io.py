import re
from collections import defaultdict
import numpy as np
import pandas as pd
from io import StringIO

from datruf_utils import (run_command,
                          add_element)

from datruf_core import (Alignment,
                         Path)


# Extract TR intervals information of reads from DBdump output
def load_dbdump(datruf):
    tr_intervals_all = defaultdict(dict)
    columns = ("dbid", "start", "end")
    counter = 0
    with open(datruf.dbdump, 'r') as f:
        for line in f:
            data = line.strip().split(' ')
            if data[0] == "R":
                dbid = int(data[1])
                if dbid > datruf.end_dbid:
                    break
            elif data[0] == "T0" and int(data[1]) > 0:
                if dbid < datruf.start_dbid:
                    continue
                for i in range(int(data[1])):
                    start = int(data[1 + 2 * (i + 1)])
                    end = int(data[1 + 2 * (i + 1) + 1])
                    add_element(tr_intervals_all, columns, counter,
                                [dbid, start, end])
                    counter += 1

    tr_intervals_all = pd.DataFrame.from_dict(tr_intervals_all)
    if len(tr_intervals_all) != 0:
        tr_intervals_all = tr_intervals_all.loc[:, columns]
    return tr_intervals_all


# Extract alignments information of reads from LAdump output
def load_ladump(datruf):
    alignments_all = defaultdict(dict)
    columns = ("dbid", "abpos", "aepos", "bbpos", "bepos")
    counter = 0
    with open(datruf.ladump, 'r') as f:
        for line in f:
            data = line.strip().split(' ')
            if data[0] == "P":
                dbid = int(data[1])
                if dbid > datruf.end_dbid:
                    break
            elif data[0] == "C":
                if dbid < datruf.start_dbid:
                    continue
                add_element(alignments_all, columns, counter,
                            [dbid] + list(map(int, data[1:5])))
                counter += 1

    alignments_all = pd.DataFrame.from_dict(alignments_all)

    alignments_all = alignments_all.loc[:, columns]
    return alignments_all


# Load TR intervals of single read
# If those of all reads are already stored, use that data
# Otherwise, call DBdump
def load_tr_intervals(datruf):
    if hasattr(datruf, "tr_intervals_all"):
        all_data = datruf.tr_intervals_all
        tr_intervals = all_data[all_data["dbid"] == datruf.read_id]
        tr_intervals = ([] if len(tr_intervals) == 0
                        else list(tr_intervals
                                  .apply(lambda x: (x["start"], x["end"]),
                                         axis=1)))
    else:
        command = ("DBdump -mtan %s %d | awk '$1 == \"T0\" {print $0}'"
                   % (datruf.db_file, datruf.read_id))
        dbdump = run_command(command).strip().split(' ')

        tr_intervals = [(int(dbdump[1 + 2 * (i + 1)]),
                         int(dbdump[1 + 2 * (i + 1) + 1]))
                        for i in range(int(dbdump[1]))]

    return tr_intervals


# Load alignments of single read like TR intervals
def load_alignments(datruf):
    if hasattr(datruf, "alignments_all"):
        all_data = datruf.alignments_all
        alignments = all_data[all_data["dbid"] == datruf.read_id]
    else:
        command = ("LAdump -c %s %s %d | awk '$1 == \"C\" {print $0}'"
                   % (datruf.db_file, datruf.las_file, datruf.read_id))
        ladump = StringIO(run_command(command))

        alignments = pd.read_csv(ladump, sep=" ",
                                 names=("abpos", "aepos", "bbpos", "bepos"))

    alignments = (alignments
                  .assign(distance=lambda x: x["abpos"] - x["bbpos"])
                  .sort_values(by="abpos", kind="mergesort")
                  .sort_values(by="distance", kind="mergesort")
                  .reset_index(drop=True))
    return alignments


def convert_symbol(aseq, bseq, symbols):
    converted = ""
    for i in range(len(symbols)):
        if symbols[i] == '|':
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
            if len(data) == 0:   # empty line due to "]"
                continue
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
        paths.append(Path(ab, ae, bb, be,
                          Alignment(aseq,
                                    bseq,
                                    convert_symbol(aseq, bseq, symbols))))

    return paths
