#!/usr/bin/env python

import sys
import argparse
import gzip
from lib import misc, tbl

"""
Read in a set of variant files (BIM, VCF, or daner). Take the 
intersect (or union) of chr-pos-varid. print varid to stdout.
"""

def main():

    args = parse_args()

    if len(args.var_files) % 2 != 0:
        print("ERROR : file format for at least one var_file not provided.")
        sys.exit(1)

    var_id_set = set()
    i = 0
    while i < len(args.var_files):
        var_file = args.var_files[i]
        var_file_format = args.var_files[i+1]

        print("Subsetting on file : " + var_file)
        print("File format : " + var_file_format)

        has_header = False
        if var_file_format in ("daner"):
            has_header = True
        
        tbl_i = tbl.Tbl(var_file, with_header=has_header)
        var_id_set_i = set()
        j = 0
        while True:
            j += 1
            var_row = tbl_i.get_row(return_dict=has_header)
            if var_row == -1:
                break
            elif var_row == None:
                continue
            
            if var_file_format == "daner":
                chr_pos_varid = (var_row["CHR"],var_row["BP"],var_row["SNP"])
            elif var_file_format == "bim":
                chr_pos_varid = (var_row[0], var_row[3], var_row[1])
            else:
                # default : look for vcf format
                chr_pos_varid = (var_row[0], var_row[1], var_row[2])
            var_id_set_i.add(chr_pos_varid)
       
        tbl_i.close_fh()

        if args.union == True:
            var_id_set = var_id_set.union(var_id_set_i)
        elif i > 0:
            var_id_set = var_id_set.intersection(var_id_set_i)
        else:
            var_id_set = var_id_set_i
            
        print("Total variants present in file : " + str(len(var_id_set_i)))
        print("Total variants in current master set : " + str(len(var_id_set)))

        i += 2

    print("Writing variant IDs to output file.")

    # init the output file 
    out_fh = open(args.varid_out_file, "w")
    out_fh.write("")
    out_fh.close()

    # write variant IDs to output file
    out_list = []
    i = 0
    for chr_pos_varid in var_id_set:
        out_list.append(chr_pos_varid[2])
        i += 1
        if i == args.write_size:
            file_append(args.varid_out_file, out_list)
            i = 0
            out_list = []
    if len(out_list) > 0:
        file_append(args.varid_out_file, out_list)

    print("Variant subsetting complete.")

    return

def file_append(filename, out_list):
    out_fh = open(filename, "a")
    out_str = "\n".join(out_list)
    out_fh.write(out_str + "\n")
    out_fh.close()
    return

def get_chr_pos_varid(var_line,
                      var_file_format):
    if var_file_format == "bim":
        data = var_line.rstrip().split()
        chr_pos_varid = (data[0], data[3], data[1])
    elif var_file_format == "vcf":
        data = var_line.rstrip().split()
        chr_pos_varid = (data[0], data[1], data[2])
    elif var_file_format == "daner":
        data = var_line.rstrip().split()
        chr_pos_varid = (data[0], data[2], data[1])
    else:
        # else, assume some vcf-styled format
        chr_pos_varid = (data[0], data[1], data[2])
    return chr_pos_varid

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--union", action="store_true", default=False,
                      help="instead of taking intersect of var ids, take union.")
    args.add_argument("--write-size", type=int,
                      action="store", default=10000, 
                      help="Number of variant IDs that go into output file in a single write.")
    args.add_argument("--tbl-delim", type=str, action="store",
                      default="\t", help="delimiter for input table file.")
    args.add_argument("--out-delim", type=str, action="store",
                      default="\t", help="delimiter for output table file.")
    args.add_argument("varid_out_file", type=str, action="store",
                      help="output file with variant IDs.")
    args.add_argument("var_files", nargs="+", 
                      help="files with variants, in format '<var_file> <file_format>. " + \
                           "allowed file formats : bim, vcf, daner.")    
    return args.parse_args()

if __name__ == "__main__":
    main()
