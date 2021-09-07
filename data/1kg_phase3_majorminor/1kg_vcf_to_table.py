#!/usr/bin/env python3

import gzip
import sys

## PARAMETERS
AMBIGUOUS_ONLY=True
AMBIGUOUS_REFALT=set(["AT","TA","GC","CG"])
INFO_FIELDS_STORE_STR="AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF"
AFMIN_COL="AF"

def main():

    # get parameters
    global AMBIGUOUS_ONLY
    global AMBIGUOUS_REFALT

    ARGS = sys.argv[1:]
    try:
        in_vcf_gz = ARGS[0]
        afmin = float(ARGS[1])
        outroot = ARGS[2]
    except:
        print("1kg_vcf_to_table.py <in.vcf.gz|stdin> <af_min> <outroot>")
        sys.exit(1)

    # get info fields
    info_fields_store_str = INFO_FIELDS_STORE_STR

    # init filehandle to vcf
    if in_vcf_gz == "stdin":
        vcf_fh = sys.stdin
    elif in_vcf_gz.find(".vcf.gz") != -1:
        vcf_fh = gzip.open(in_vcf_gz, "rt")
    else:
        vcf_fh = open(in_vcf, "r")

    # get info fields
    info_fields_store = info_fields_store_str.split(",")
    info_fields_store_set = set(info_fields_store)

    # print header to stdout
    header_list = ["CHROM","POS","REF","ALT"] + info_fields_store
    header_str = "\t".join(header_list)

    # define previous chrom
    prev_chrom=None

    # init array of lines to write to chromosome-level files
    out_lines = [header_str]

    # for each line..
    for line in vcf_fh:

        # skip if commented out 
        if line[0] == "#": continue

        # get data from line
        data = line.rstrip().split("\t")
        chrom = data[0].replace("chr","")
        pos = int(data[1])
        ref = data[3]
        alts_str = data[4]

        # if on a new chromosome, then time to write a new file
        if chrom != prev_chrom and prev_chrom != None:
            out_tsv_gz = outroot + ".chr" + prev_chrom + ".tsv.gz"
            write_tsv_gz(out_lines, out_tsv_gz)
            out_lines=[header_str]

        # update prev chromosome
        prev_chrom = chrom

        # get info
        info_str = data[7]
        info_list = info_str.split(";")
        info_dict = dict()
        for info_keyval_str in info_list:
            info_keyval = info_keyval_str.split("=")
            if len(info_keyval) != 2: continue
            info_key = info_keyval[0]
            info_vals = info_keyval[1]
            if info_key in info_fields_store_set:
                info_dict[info_key] = info_vals.split(",")

        # for each alt ..
        alts = alts_str.split(",")
        for i in range(len(alts)):
            alt = alts[i]


            # skip if focusing on ambiguous only and var is nonambiguous
            if AMBIGUOUS_ONLY == True:
                if ref+alt not in AMBIGUOUS_REFALT:
                    continue

            # skip if not found at afmin
            afmin_val = float(info_dict[AFMIN_COL][i])
            if afmin_val <= afmin or afmin_val >= (1-afmin): continue

            vals = []
            for info_key in info_fields_store:
                if info_key in info_dict:
                    af_val = float(info_dict[info_key][i])
                    if af_val < 0.5:
                        val = "2"
                    elif af_val > 0.5:
                        val = "1"
                    else:
                        val = "0"
                else:
                    val = "2"
                vals.append(str(val))

            out_row_list=[chrom, str(pos), ref, alt] + vals
            out_row = "\t".join(out_row_list)
            out_lines.append(out_row)

    # one last file write
    out_tsv_gz = outroot + ".chr" + chrom + ".tsv.gz"
    write_tsv_gz(out_lines, out_tsv_gz)

    return
 
def write_tsv_gz(out_lines, out_filename):
    print("WRITING : " + out_filename)
    out_fh = gzip.open(out_filename, "wt")
    for line in out_lines:
        out_fh.write(line + "\n")
    out_fh.close()
    return
    
if __name__ == "__main__":
    main()
