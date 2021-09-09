
import sys
import gzip
import io
import math

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
        fh_gz = gzip.open(filename, "rb")
        fh = io.BufferedReader(fh_gz)
    else:
        fh = open(filename, "r")
    return fh

def load_fam_trios(fam_filename):
    proband_parents = {}
    fh = open(fam_filename,"r")
    for line in fh:
        data = line.rstrip().split()
        if data[5] == "2":
            proband_parents[data[1]] = [data[2],
                                        data[3]]
    fh.close()
    return proband_parents

def between(val_x,range_y):
    if range_y[0] <= val_x and val_x <= range_y[1]:
        return True
    else:
        return False

def beta_to_or(beta_val):
    try:
        x = math.exp(float(beta_val))
    except:
        x = "NA"
    return x

def or_to_beta(or_val):
    try:
        x = math.log(float(or_val))
    except:
        x = "NA"
    return x

def ref_fh_is_withchr(ref_fh, ref_fh_type):
    # assumed that generally most ref files will be without chr
    is_withchr=False
    if ref_fh_type == "pysam":
        chroms = ref_fh.references
        for chrom in chroms:
            if chrom.find("chr") == 0:
                is_withchr = True
    elif ref_fh_type == "twobitreader":
        chroms =  ref_fh.keys()
        for chrom in chroms:
            if chrom.find("chr") == 0:
                is_withchr = True
    else:
        raise ValueError("filehandle type '"+ref_fh_type+"' not supported.")
    return is_withchr
