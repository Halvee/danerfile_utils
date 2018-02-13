
#                    if args.max_csq == True:                                    
#                    for max_col in MAX_COLS:                                
#                    full_header.extend(["CSQ_MAX_"+max_col])

import sys
import argparse
from genetics_munge_lib import tbl,misc,annot,vcf

def main():
    args = parse_args()
    csq_names = args.csq_names.split(",")
    tbl_x = tbl.Tbl(args.tbl_file, delim=args.tbl_delim)
    i = 0
    while(1):
        i += 1
        tbl_x.get_row(return_dict=True)
        if i == 1: 
            for csq_str_i in args.csq_threshes:
                csq_list_i = csq_str_i.split(":") 
                csq_names_out_i = csq_list_i[0].split(",")
                tbl_x.header_list.extend(csq_names_out_i)
            print(args.tbl_delim.join(tbl_x.header_list))
            continue
        
        if len(tbl_x.row_list) == 0 or len(tbl_x.row_dict) == 0: break
        csq=annot.AnnotTxs(csq_names,  
                           tbl_x.row_dict[args.csq_colname])
        for csq_str_i in args.csq_threshes:
            csq_list_i = csq_str_i.split(":")
            csq_cols_return = csq_list_i[1].split(",")
            csq_key = csq_list_i[2]
            csq_func = csq_list_i[3]
            out = csq.max_csq(csq_cols_return, csq_key, csq_func)
            tbl_x.row_list.extend(out)

        print("\t".join(tbl_x.row_list))
    return    

def get_max_csqs(dat_row, vcf_r):                                               
    in_str=vcf_r.metainfo_dicts["CSQ"]["Description"]                           
    csq_header=annot.process_csq_header_desc(in_str)                            
    csq=annot.AnnotTxs(csq_header,                                              
                       vcf_r.vcf_entry.info_dict["CSQ"])                        
    max_csq = csq.max_csq(sift_min=True,polyphen_max=True)                      
    for max_col in MAX_COLS:                                                    
        dat_row.extend([max_csq["CSQ_MAX_"+max_col]])                           
    return dat_row

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--tbl-delim", type=str,
                      action="store", default="\t", 
                      help="Delimiter for seperating columns.")
    args.add_argument("--cols-replace-info",type=str,
                      action="store",default=None,
                      help="comma delimited set of old:new column names for INFO values.")
    args.add_argument("--cols-replace-format",type=str,
                      action="store",default=None,
                      help="comma delimited set of old:new column " + \
                           "names for FORMAT values.")
    args.add_argument("--csq-colname", type=str, 
                      action="store", default="CSQ",
                      help="colname for consequence data in table.")
    args.add_argument("tbl_file", type=str,
                      help="input table file.")
    args.add_argument("csq_names", type=str,
                      help="comma delim names of CSQ items")
    args.add_argument("csq_threshes", nargs="+",
                      help="csq_names_report:csq_names_test:csq_name_test:max/min, " + \
                           "for each csq_name_report you want a col for.")
    return args.parse_args()

if __name__ == "__main__":
    main()
