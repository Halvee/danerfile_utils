#!/usr/bin/env python

import sys
import argparse
import gzip
from danerfile_utils_lib import misc, tbl

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
            tbl_i.get_row(return_dict=has_header)
            if tbl_i.row_str == "":
                break
            elif tbl_i.row_str[0] == "#":
                continue
            
            if var_file_format == "daner":
                id = tbl_i.row_dict["SNP"]
                if args.chr_pos_varid == True:
                    id = (tbl_i.row_dict["CHR"],
                          tbl_i.row_dict["BP"],
                          tbl_i.row_dict["SNP"])
            elif var_file_format == "bim":
                id = tbl_i.row_list[1]
                if args.chr_pos_varid == True:
                    id = (tbl_i.row_list[0], tbl_i.row_list[3], tbl_i.row_list[1])
            elif var_file_format == "list":
                if args.chr_pos_varid == True:
                    print("ERROR : 'list' file format cannot contain chrom or pos columns.")
                    sys.exit(1)
                id = tbl_i.row_list[0]
            else:
                # default : in vcf fashion, expect format chr, bp, var_id
                id = tbl_i.row_list[2]
                if args.chr_pos_varid == True:
                    id = (tbl_i.row_list[0], tbl_i.row_list[1], tbl_i.row_list[2])
            if id in var_id_set or i == 0:
                var_id_set_i.add(id)
       
        tbl_i.close_fh()

        var_id_set = var_id_set_i
            
        print("Total variants present in file : " + str(j))
        print("Total intersecting variants : " + str(len(var_id_set)))

        i += 2

    print("Writing variant IDs to output file.")

    # init the output file 
    out_fh = open(args.varid_out_file, "w")
    out_fh.write("")
    out_fh.close()

    # write variant IDs to output file
    out_list = []
    i = 0
    for id in var_id_set:
        if args.chr_pos_varid == True:
            id_out = id[2]
        else:
            id_out = id
        out_list.append(id_out)
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

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--write-size", type=int,
                      action="store", default=10000, 
                      help="Number of variant IDs that go into output file in a single write.")
    args.add_argument("--tbl-delim", type=str, action="store",
                      default="\t", help="delimiter for input table file.")
    args.add_argument("--out-delim", type=str, action="store",
                      default="\t", help="delimiter for output table file.")
    args.add_argument("--chr-pos-varid", action="store_true",default=False,
                      help="store chr,pos,varid as key in comparing variant IDs.")
    args.add_argument("varid_out_file", type=str, action="store",
                      help="output file with variant IDs.")
    args.add_argument("var_files", nargs="+", 
                      help="files with variants, in format '<var_file> <file_format>. " + \
                           "allowed file formats : bim, vcf, daner, list.")    
    return args.parse_args()

if __name__ == "__main__":
    main()
