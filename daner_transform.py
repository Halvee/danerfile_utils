import argparse
from lib.tbl import Tbl, Cnds

def main():
    args = parse_args()

    """
    Init instance of Cnds class, reading in cnds file if defined
    """
    cnds = None
    if args.cnds_file != None:
        cnds = Cnds(args.cnds_file)

    """
    create instance of Tbl class for input daner file. 
    """
    daner = Tbl(args.daner_file, with_header=True, 
                cols_recode_str=args.header_cols_recode)

    """
    get and print daner header
    """
    if args.header_subset != None:
        header_list = args.header_subset.split(",")
    else:
        header_list = daner.header_list

    """
    if header recodings are defined by user, perform col name recoding
    """
    if args.header_cols_recode != None:
        header_cols_recode = args.header_cols_recode.split(",")
        col_oldnew_dict = {}
        for col_oldnew in header_cols_recode:
            (col_old, col_new) = col_oldnew.split(":")[:2]
            col_oldnew_dict[col_old] = col_new
        for i in range(len(header_list)):
            if header_list[i] in col_oldnew_dict:
                header_list[i] = col_oldnew_dict[header_list[i]]

    """
    if OR / BETA conversion occuring, switch header values
    """
    for i in range(len(header_list)):
        if args.or_to_beta == True and header_list[i] == "OR":
            header_list[i] = "BETA"
        elif args.beta_to_or == True and header_list[i] == "BETA":
            header_list[i] = "OR"

    header_str = args.delim.join(header_list)
    print(header_str)
    
    """
    iterate through daner rows, only print those that pass cnds
    """
    while(1):
        row_dict = daner.get_row(delim=None, 
                                 return_dict=True)
        if len(row_dict) == 0: 
            break
        cnds_pass = True
        if cnds != None:
            cnds_pass = cnds.test(row_dict)
        if cnds_pass == True:
            row_list = []
            
            """
            reconstruct row list
            """
            for col in header_list:
                row_list.append(str(row_dict[col]))
            row_str = args.delim.join(row_list)
            print(row_str)
            
    daner.close_fh()

    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("-d","--delim", type=str,
                      action="store", default="\t", 
                      help="Delimiter for output daner.")
    args.add_argument("--or-to-beta", action="store_true", default=False,
                      help="convert odds ratio values and corresponding std errors to beta values.")
    args.add_argument("--beta-to-or", action="store_true", default=False,
                      help="convert beta values and corresponding std errors to odds ratio values.")
    args.add_argument("--cnds-file", type=str, action="store", default=None,
                      help="File where each row is a condition that each daner row must meet, "+ \
                           "each row has format 'parameter operator threshold'")
    args.add_argument("--header-cols-recode", type=str, action="store",
                      default=None, help="comma-delimted list of old_col_name:new_col_name strings")
    args.add_argument("--header-subset", type=str, action="store",
                      default=None, help="comma-delimted list of header columns to subset on")
    args.add_argument("--dummy-cols", type=str, action="store",
                      default=None, help="comma-delimited list of dummy header columns to add, " + \
                                         "of form header_col:global_col_value")
    args.add_argument("daner_file", 
                      help="Source LD clumped file for case/control, made by PLINK, " + \
                           "with min(P) variant tagged.")
    
    return args.parse_args()

if __name__ == "__main__":
    main()

