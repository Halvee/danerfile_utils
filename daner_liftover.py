#!/usr/bin/env python

import os
import sys
import argparse
import re
import subprocess
from danerfile_utils_lib.tbl import Tbl

def main():

    # parse args
    args = parse_args()

    # check to make sure that input_build -> output_build liftover is 
    # currently supported based on the presence of obtainable chain files
    liftover_chainfile = get_liftover_chainfile(args.input_build, 
                                                args.output_build, 
                                                args.liftover_chain_dir)
    

    """
    create instance of Tbl class for input daner file. 
    """
    daner = Tbl(args.daner_file,
                delim=args.in_delim,
                with_header=True)
    
    # form bed temp input and output files, along with temp unmapped file
    bed_in = args.liftover_tempfile_dir + "/" + \
             "input." + args.input_build + ".bed"
    bed_out = args.liftover_tempfile_dir + "/" + \
              "output." + args.output_build + ".bed"   
    unmapped_txt = args.liftover_tempfile_dir + "/" + \
                   "output." + args.output_build + ".unmapped.txt"
    out_fh = open(bed_in, "w")

    """
    iterate through daner rows, and for each one, 
    get info needed for remapping
    """
    header_list = daner.header_list
    header_str = args.out_delim.join(header_list)
    while(1):
        daner.get_row(return_dict=True)
        if daner.row_str == "":
            break
        elif daner.row_str[0] == "#":
            continue
        chrom_i = daner.row_dict[args.chrom_col]
        pos_i = daner.row_dict[args.pos_col]
        snpid_i = daner.row_dict[args.snpid_col]
        end_i = int(pos_i)
        start_i = end_i - 1

        # make sure chromosome has 'chr' in front
        if chrom_i.find("chr") != 0:
            chrom_i = "chr" + chrom_i

        out_str = "\t".join([chrom_i, str(start_i), str(end_i), snpid_i]) + "\n"
        out_fh.write(out_str)
    out_fh.close()

    # execute liftover command
    sys.stderr.write("liftover exe path : " + args.liftover_exe_path + "\n")
    sys.stderr.write("input bed path : " + bed_in + "\n")
    sys.stderr.write("liftover chainfile : " + liftover_chainfile + "\n")
    sys.stderr.write("output bed path : " + bed_out + "\n")
    sys.stderr.write("putput unmapped path : " + unmapped_txt + "\n")
    subprocess.call([args.liftover_exe_path,
                     bed_in,
                     liftover_chainfile,
                     bed_out,
                     unmapped_txt])

    # read lifted-over coordinates into dictionary of 
    # snpid -> (chrom_liftedover, pos_liftedover)
    chrpos_liftover = dict()
    in_fh = open(bed_out, "r")
    for line in in_fh:
        data = line.rstrip().split()
        chr_i = data[0]
        pos_i = data[2]
        snpid_i = data[3]

        # strip 'chr' string from front of chromosome id
        chr_i = chr_i.replace("chr","")

        chrpos_liftover[snpid_i] = (chr_i, pos_i)
    in_fh.close()

    # init filehandle to output file to write 
    if args.output_daner_file == "stdout":
        out_fh = sys.stdout
    elif re.search(".gz$", args.output_daner_file, flags=0) != None:
        out_fh = gzip.open(args.output_daner_file, "wt")
    else:
        out_fh = open(args.output_daner_file, "w")

    # set byte index in original danerfile to line 1 of data
    daner.fh.seek(0)
    line1 = daner.fh.readline()

    # reread input danerfile, remapping chrom/pos to liftered-over coordinates
    i0=0
    i=0
    header_list = daner.header_list
    header_str = args.out_delim.join(header_list)
    out_fh.write(header_str + "\n")
    while(1):
        i0+=1
        daner.get_row(return_dict=True)
        if daner.row_str == "":
            break
        elif daner.row_str[0] == "#":
            continue
        chrom_i = daner.row_dict[args.chrom_col]
        pos_i = daner.row_dict[args.pos_col]
        snpid_i = daner.row_dict[args.snpid_col]
        end_i = int(pos_i)
        start_i = end_i - 1
        
        # skip entry if snpid not successfully lifted over
        if snpid_i not in chrpos_liftover: continue

        # get lifted-over chrom, pos for snpid
        chrom_liftover_i = chrpos_liftover[snpid_i][0]
        pos_liftover_i = chrpos_liftover[snpid_i][1]

        # overwrite the original chrom, pos entries in snpfile
        daner.row_dict[args.chrom_col] = chrom_liftover_i
        daner.row_dict[args.pos_col] = pos_liftover_i


        """
        reconstruct row list
        """
        row_list = []
        for col in header_list:
            row_list.append(str(daner.row_dict[col]))
        row_str = args.out_delim.join(row_list)

        # write to output filehandle
        out_fh.write(row_str + "\n")
        i += 1

    # close filehandle to new snpfile, original snpfile
    out_fh.close()
    daner.fh.close()

    # unless indicated otherwise by user delete temporary files
    if args.keep_liftover_tempfiles == False:
        for x in (bed_in, bed_out, unmapped_txt):
            subprocess.call(['rm', x])

    return

def get_liftover_chainfile(build_from, build_to, liftover_chain_dir):
    liftover_chainfile=liftover_chain_dir + \
                       build_from + "To" + \
                       build_to.replace("h","H") + ".over.chain.gz"
    if os.path.isfile(liftover_chainfile) == False:
        print("ERROR : ")
        print(build_from + " -> " + build_to + " liftover not supported. Need corresponding chain file.")
        dir_contents = os.listdir(liftover_chain_dir)  
        dir_contents.sort(reverse=True)
        print("Available conversion supported based on chain files :")
        for item in dir_contents:
            if item.find(".over.chain.gz") != -1:
                item = item.replace(".over.chain.gz","")
                item = item.replace("To"," -> ")
                item = item.replace("H","h")
                print(item)
        sys.exit(1)
    return liftover_chainfile

def chainfile_absence_msg(build_from_to):
    global HG_CHAINFILES
    print("ERROR : ")
    print(build_from_to + " not supported.")
    print("The following liftOver operations are supported currently :")
    for hg_from_to in HG_CHAINFILES:
        print(hg_from_to)
    sys.exit(1)
    return

def parse_args():
    parser = argparse.ArgumentParser(description='Use UCSC genome liftover executable to convert daner coords.')
    parser.add_argument('--liftover-exe-path', type=str, action='store',
                        default=sys.path[0] + '/bin/liftOver',
                        help='path to ucsc genome liftOver utility')
    parser.add_argument('--liftover-chain-dir', action='store', type=str,
                        default=sys.path[0]+'/liftOver_chain_files/',
                        help='path to directory with liftOver chain files')
    parser.add_argument('--liftover-tempfile-dir', action='store', type=str,
                        default='/tmp/',
                        help='path to directory where temporary BED ' + \
                             'files will be written')
    parser.add_argument('--keep-liftover-tempfiles', action='store_true',
                        default=False,
                        help='keep liftover tempfiles, rather than ' + \
                             'removing them by default.')
    parser.add_argument('--snpid-col', type=str, default='SNP',
                        help='name of column with SNP ID.')
    parser.add_argument('--chrom-col', type=str, default='CHR',
                        help='name of column with chromosome.')
    parser.add_argument('--pos-col', type=str, default='BP',
                        help='name of column with position')
    parser.add_argument('input_build', type=str,
                        help='reference genome build for coords in input file')
    parser.add_argument('--in-delim', type=str, default="\t",
                        help='Delimiter for input daner.')
    parser.add_argument('--out-delim', type=str, default="\t",
                        help='Delimiter for output daner.')
    parser.add_argument('daner_file', type=str,
                        help="Input daner file. Write 'stdin' to read ' + \
                             'from standard input")
    parser.add_argument('output_build', type=str,
                        help='reference genome guild for coords in output file')
    parser.add_argument('output_daner_file', type=str,
                        help="output daner file with repmapped " + \
                             "coordinates. Write 'stdout' to write " + \
                             "to standard output.")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()
