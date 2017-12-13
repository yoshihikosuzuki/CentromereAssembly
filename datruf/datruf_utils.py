import sys
import subprocess
from interval import interval


def run_command(command):
    try:
        out = subprocess.check_output(command, shell=True)#,
                                      #stderr=subprocess.STDOUT)   # TODO: fix LAdump's bug and then comment-in
    except subprocess.CalledProcessError as proc:
        print(proc.output.decode('utf-8'))
        sys.exit(1)
    else:
        return out.decode('utf-8')


# Add an element (corresponding to a row in pd.DataFrame)
# into a defaultdict(dict) that will be in the end converted to pd.DataFrame
def add_element(dic, columns, counter, data):
    for i in range(len(columns)):
        dic[columns[i]][counter] = data[i]


# Generate a line shape object in plotly
def make_line(x0, y0, x1, y1, col, width):
    return {'type': 'line', 'xref': 'x', 'yref': 'y',
            'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1,
            'line': {'color': col, 'width': width}}


def interval_len(intvls):
    ret = 0
    for intvl in intvls.components:
        ret += intvl[0][1] - intvl[0][0] + 1
    return ret


# subtraction of integer intervals: A - B, and afterwards remove intervals
# whose length is less than "length_threshold"
def subtract_interval(a_interval, b_interval, length_threshold=0):
    ret_interval = interval()
    intersection = a_interval & b_interval
    a_index = 0
    i_index = 0
    # use when there might be other intersecting intervals in an A-interval
    flag_load_new_interval = True
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
                    # the tail interval that did not intersect
                    # with this B-interval (but maybe with the next B-interval)
                    a_intvl = interval[i_end + 1, a_end]
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
