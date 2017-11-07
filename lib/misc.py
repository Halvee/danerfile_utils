
import sys
import gzip

def keyval_list_pair_to_dict(list_i, list_j):
    y = {}
    if len(list_i) != len(list_j):
        raise Exception("list_i and list_j are unequal length.")
    for x in range(len(list_i)):
        y_key = list_i[x]
        y_val = list_j[x]
        y[y_key] = y_val
    return y

def open_file(filename):
    if filename == "stdin":
        fh = sys.stdin
    elif filename.find(".gz") != -1:
        fh = gzip.open(filename, "rb")
    else:
        fh = open(filename, "r")
    return fh
