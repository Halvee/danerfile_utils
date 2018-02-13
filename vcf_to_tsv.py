
import argparse
import io
from genetics_munge_lib import vcf,tbl,misc,annot

## GLOBAL VARS
MAX_COLS=["IMPACT","IMPACT_GENESYMBOL","SIFT","PolyPhen"]

def main():
    args = parse_args()

    filter_classifs = None
    if args.filter_classifs != None:
        filter_classifs = set(args.filter_classifs.split(","))

    gts_exclude = None
    if args.gts_exclude != None:
        gts_exclude = set(args.gts_exclude.split(","))

    printed_header = False

    for vcf_file in args.vcf_files:
        
        vcf_fh = misc.open_file(vcf_file)
        vcf_f = io.BufferedReader(vcf_fh)

        vcf_r = vcf.VcfReader(vcf_f, cols_replace_format=args.cols_replace_format,
                              cols_replace_info=args.cols_replace_info)

        while(1):
            dat = vcf_r.next_line()
            if dat[1] == 0:
                pass
            elif dat[1] == 1:
                if printed_header == False:
                    
                    full_header = ["VARIANT_ID", "SAMPLE_ID"] + vcf_r.vcf_header[:7] + \
                                  vcf_r.metainfo_lists["INFO"] + \
                                  vcf_r.metainfo_lists["FORMAT"]
                    print(args.delim.join(full_header))
                    printed_header = True
            elif dat[1] == 2:
                x = vcf_r.vcf_entry
                if filter_classifs != None: 
                    if x.filter not in filter_classifs: continue
                dat_rows = x.get_sample_rows(vcf_r.metainfo_lists,
                                             vcf_r.sample_list,
                                             gts_exclude=gts_exclude)
                for dat_row in dat_rows:
                    dat_row = [str(val) for val in dat_row]
                    print("\t".join(dat_row))
            else:
                break

    return

def cols_replace(cols_list, cols_replace):
    if cols_replace==None:
        return cols_list
    cols_replace = cols_replace.split(",")
    cols_replace = [i.split(":") for i in cols_replace]
    cols_replace_d = {i[0]:i[1] for i in cols_replace}
    for i in range(len(cols_list)):
        if cols_list[i] in cols_replace_d:
            old = cols_list[i]
            cols_list[i] = cols_replace_d[old]

    return cols_list


def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--delim", type=str,
                      action="store", default="\t", 
                      help="Delimiter for seperating columns.")
    args.add_argument("--cols-replace-info",type=str,
                      action="store",default=None,
                      help="comma delimited set of old:new column names for INFO values.")
    args.add_argument("--cols-replace-format",type=str,
                      action="store",default=None,
                      help="comma delimited set of old:new column " + \
                           "names for FORMAT values.")
    args.add_argument("--filter-classifs",type=str,
                      action="store",default=None)
    args.add_argument("--gts-exclude", type=str, default=None, action="store",
                      help="genotypes to exclude from output tsv.")
    args.add_argument("vcf_files", nargs="+",
                      help="input VCF files.")
    return args.parse_args()

if __name__ == "__main__":
    main()
