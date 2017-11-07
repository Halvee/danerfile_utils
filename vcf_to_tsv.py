
import argparse
from lib import vcf,tbl,misc

def main():
    args = parse_args()
    
    printed_header = False

    cnds = tbl.Cnds(args.cnds_file)

    for vcf_file in args.vcf_files:
        
        vcf_fh = misc.open_file(vcf_file)
        
        vcf_r = vcf.VcfReader(vcf_fh)

        while(1):
            dat = vcf_r.next_line()
            if dat[1] == 0:
                pass
            elif dat[1] == 1:
                if printed_header == False:
                    full_header = vcf_r.vcf_header[:7] + \
                                  vcf_r.metainfo_lists["INFO"] + \
                                  vcf_r.metainfo_lists["FORMAT"]
                    print(args.delim.join(full_header))
                    printed_header = True
            elif dat[1] == 2:
                dat_rows = vcf_r.vcf_entry.get_sample_rows(vcf_r.metainfo_lists, vcf_r.sample_list)
                for dat_row in dat_rows:
                    if len(cnds.cnds) > 0:
                        dat_dict = misc.keyval_list_pair_to_dict(full_header,dat_row)
                        if cnds.test(dat_dict) == True:
                            dat_row = [str(val) for val in dat_row]
                            print("\t".join(dat_row))
            else:
                break

    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--delim", type=str,
                      action="store", default="\t", 
                      help="Delimiter for seperating columns.")
    args.add_argument("--eff-type", type=str, choices=("OR","BETA"),
                      action="store", default="OR",
                      help="Effect direction type to use, corresponds to col names.")
    args.add_argument("--cnds-file", type=str, action="store",
                      default=None, help="file with conditions to subset on.")
    args.add_argument("vcf_files", nargs="+",
                      help="input VCF files.")
    
    return args.parse_args()

if __name__ == "__main__":
    main()
