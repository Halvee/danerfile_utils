
import argparse
from lib import tbl,misc

def main():
    args = parse_args()
    
    tbl_i = tbl.Tbl(args.tbl_file, delim=args.tbl_delim)

    if len(args.cnds_in_out_files) % 2 != 0:
        print("ERROR : need an output file name for each input cnds file")
        sys.exit(1)

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
        row_dict, row_list = tbl_i.get_row(return_dict=True,
                                           return_list_too=True)
        if len(row_list) == 0 or len(row_dict) == 0: break
        for i in range(len(cnds_sets)):
            if cnds_sets[i].test(row_dict) == True:
                row_str = args.out_delim.join(row_list)
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
    args.add_argument("tbl_file", type=str, action="store",
                      help="name of input table file.")
    args.add_argument("cnds_in_out_files", nargs="+",
                      help="cnds files and corresponding output files (cnds_1, out_1, cnds_2, out_2...).")
    
    return args.parse_args()

if __name__ == "__main__":
    main()
