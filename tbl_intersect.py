
import os
import sys
import argparse
from danerfile_utils_lib import tbl

def main():
    args = parse_args()
    
    intersect_set = set()
    fh = open(args.intersect_file, "r")
    for line in fh: intersect_set.add(line.rstrip())
    fh.close()

    tbl_i = tbl.Tbl(args.tbl_file, delim=args.delim, with_header=args.has_header)

    if args.has_header: print(tbl_i.header_str)

    while(1):
        tbl_i.get_row(return_dict=args.has_header)
        if tbl_i.row_str == "": break
        try:
            key_i = tbl_i.row_dict[args.tbl_col]
        except:
            key_i = tbl_i.row_list[int(args.tbl_col)]
        
        if key_i not in intersect_set: continue
        print(tbl_i.row_str.rstrip())

    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--delim", type=str,
                      action="store", default=None, 
                      help="Delimiter for seperating columns.")
    args.add_argument("--has-header", action="store_true",default=False,
                      help="input table has header, take into account " + \
                           "when reading and writing to stdout.")
    args.add_argument("tbl_file", type=str, action="store",
                      help="name of input table file.")
    args.add_argument("tbl_col", type=str, action="store",
                      help="column name or number to intersect on.")
    args.add_argument("intersect_file", type=str, action="store",
                      help="list of values to intersect tbl_file on, at tbl_col.")
    
    return args.parse_args()

if __name__ == "__main__":
    main()
