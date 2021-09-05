
import os
import sys
import argparse
from danerfile_utils_lib import tbl,misc,annot,vcf

def main():
    args = parse_args()
    
    tbl_i = tbl.Tbl(args.tbl_file, delim=args.tbl_delim)

    if len(args.cnds_in_out_files) % 2 != 0:
        print("ERROR : need an output file name for each input cnds file")
        sys.exit(1)

    if args.min_perc_alt > 0.0:
        tbl_i.header_list.append("PERC_ALT")    
    if args.max_impact != None:
        tbl_i.header_list.extend(["GENE_NAME","IMPACT"])

    i = 0 
    cnds_sets = []
    out_files = []
    while i < len(args.cnds_in_out_files):
        cnds_file = args.cnds_in_out_files[i]
        out_file = args.cnds_in_out_files[i+1]
        cnds = tbl.Cnds(cnds_file)
        cnds_sets.append(cnds)
        fh = open(out_file, "w")
        fh.write(args.out_delim.join(tbl_i.header_list) + "\n")
        fh.close()
        out_files.append(out_file)
        i += 2
    while(1):
        tbl_i.get_row(return_dict=True)
        if len(tbl_i.row_list) == 0 or len(tbl_i.row_dict) == 0: break
        if args.min_perc_alt > 0.0:
            tbl_i.row_dict["PERC_ALT"] = vcf.ad_min_perc_alt(row_dict["AD"])
            if tbl_i.row_dict["PERC_ALT"] < args.min_perc_alt: continue
            tbl_i.row_dict["PERC_ALT"] = str(tbl_i.row_dict["PERC_ALT"])
            tbl_i.row_list.append(row_dict["PERC_ALT"])
    
        if args.max_impact != None:
            annot_txs = annot.SnpeffAnnotTxs(row_dict[args.max_impact],
                                             format=args.max_impact)
            annot_txs.get_max_impact(protein_coding_only = args.protein_coding_only)
            if annot_txs.max_eff == None: continue
            tbl_i.row_dict["GENE_NAME"] = annot_txs.max_eff.gene_name
            tbl_i.row_dict["IMPACT"] = annot_txs.max_eff.impact
            tbl_i.row_list.extend([tbl_i.row_dict["GENE_NAME"], 
                                   tbl_i.row_dict["IMPACT"]])
        for i in range(len(cnds_sets)):
            if cnds_sets[i].test(tbl_i.row_dict) == True:
                row_str = args.out_delim.join(tbl_i.row_list)
                fh = open(out_files[i], "a")
                fh.write(row_str + "\n")
                fh.close()

    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--delim", type=str,
                      action="store", default="\t", 
                      help="Delimiter for seperating columns.")
    args.add_argument("--tbl-delim", type=str, action="store",
                      default="\t", help="delimiter for input table file.")
    args.add_argument("--out-delim", type=str, action="store",
                      default="\t", help="delimiter for output table file.")
    args.add_argument("--min-perc-alt", type=float, action="store",
                      default=0.0, help="min percent of reads that are alt allele.")
    args.add_argument("--max-impact", type=str, action="store", choices=["ANN"],
                      default=None, help="convert ANN or EFF to max impact, add GENE_NAME and IMPACT")
    args.add_argument("--protein-coding-only", action="store_true", default=False,
                      help="consider annotations from protein coding transcripts only.")
    args.add_argument("tbl_file", type=str, action="store",
                      help="name of input table file.")
    args.add_argument("cnds_in_out_files", nargs="+",
                      help="cnds files and corresponding output files (cnds_1, out_1, cnds_2, out_2...).")
    
    return args.parse_args()

if __name__ == "__main__":
    main()
