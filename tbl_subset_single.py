
import os
import sys
import argparse
from genetics_munge_lib import tbl,misc,annot,vcf

def main():
    args = parse_args()
    
    tbl_i = tbl.Tbl(args.tbl_file, delim=args.tbl_delim)

    cnds = tbl.Cnds(args.cnds_file)

    i = 0
    while(1):
        i += 1
        tbl_i.get_row(return_dict=True)
        if i == 1:
            row_str = args.tbl_delim.join(tbl_i.header_list)
            print(row_str)
            continue
        if len(tbl_i.row_list) == 0 or len(tbl_i.row_dict) == 0: break
        if args.min_perc_alt > 0.0:
            tbl_i.row_dict["PERC_ALT"] = vcf.ad_min_perc_alt(row_dict["AD"])
            if tbl_i.row_dict["PERC_ALT"] < args.min_perc_alt: continue
            tbl_i.row_dict["PERC_ALT"] = str(tbl_i.row_dict["PERC_ALT"])
            tbl_i.row_list.append(row_dict["PERC_ALT"])
    
        if cnds.test(tbl_i.row_dict) == True:
            row_str = args.tbl_delim.join(tbl_i.row_list)
            print(row_str)

    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--delim", type=str,
                      action="store", default="\t", 
                      help="Delimiter for seperating columns.")
    args.add_argument("--tbl-delim", type=str, action="store",
                      default="\t", help="delimiter for input table file.")
    args.add_argument("--min-perc-alt", type=float, action="store",
                      default=0.0, help="min percent of reads that are alt allele.")
    args.add_argument("tbl_file", type=str, action="store",
                      help="name of input table file.")
    args.add_argument("cnds_file", type=str, action="store",
                      help="input cnds file")
    
    return args.parse_args()

if __name__ == "__main__":
    main()
