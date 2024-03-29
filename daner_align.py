import argparse
import sys
import re
from danerfile_utils_lib.tbl import Tbl, Cnds
from danerfile_utils_lib.daner import Marker
from danerfile_utils_lib.misc import ref_fh_is_withchr
from danerfile_utils_lib.poprefalt import PopRefAlt, a1_a2_direction

def main():
    args = parse_args()

    """
    Load reference file data, determine if contig names have 'chr' in front
    """
    searchres_fasta_gz=re.search(".fasta.gz$", args.ref_file, flags=0)
    searchres_fa_gz=re.search(".fa.gz$", args.ref_file, flags=0)
    searchres_fasta=re.search(".fasta$", args.ref_file, flags=0)
    searchres_fa=re.search(".fa$", args.ref_file, flags=0)
    searchres_2bit=re.search(".2bit$", args.ref_file, flags=0)
    is_2bit=False
    if searchres_fasta_gz != None or searchres_fa_gz != None or \
       searchres_fasta != None or searchres_fa != None:
        try:
            import pysam
        except:
            sys.exit("\nERROR : module 'pysam' required to parse ref fasta file.\n")
        ref_fh = pysam.FastaFile(args.ref_file)
        is_withchr = ref_fh_is_withchr(ref_fh, "pysam")
    elif searchres_2bit != None:
        try:
            import twobitreader
        except:
            sys.exit("\nERROR : module 'twobitreader' required to parse ref 2bit file.\n")
        ref_fh = twobitreader.TwoBitFile(args.ref_file)
        is_2bit=True
        is_withchr = ref_fh_is_withchr(ref_fh, "twobitreader")
    else:
        sys.exit("\nERROR : ref_file needs to be one of the following formats : \n")

    """
    create instance of Tbl class for input daner file. 
    """
    daner = Tbl(args.daner_file,
                delim=args.in_delim,
                with_header=True)
    if len(daner.header_list) <= 1:
        sys.exit("\nERROR : Only 0 or 1 columns read in. Do you have --in-delim set right?\n")

    """
    Use BETA as effsize column by default. If absent, test for presence of OR.
    If OR column is absent, raise exception.
    """
    effsize_col = "BETA"
    if "BETA" not in daner.header_list:
        effsize_col = "OR"
        if "OR" not in daner.header_list:
            sys.exit("\nERROR : required effsize column (BETA or OR) not found in danerfile.\n")

    """
    if defined, get formatted ambiguous freq range where vars inside are excluded
    """
    ambiguous_af_range_exclude = None
    if args.ambiguous_af_range_exclude != None:
        excl_list = args.ambiguous_af_range_exclude.split("-")
        try:
            ambiguous_af_range_exclude = [float(excl_list[0]),float(excl_list[1])]
        except:
            sys.exit("ERROR : input for opt --ambiguous-af-range-exclude not " \
                     "in format 'LOWERLIMIT-UPPERLIMIT'")

    """
    if defined, load pop refalt tsv
    """
    pop_refalt = None
    if args.pop_refalt_prefix != None:
        # init refalt balance obj
        pop_refalt = PopRefAlt(
                               delim="\t",
                               pop_refalt_af_col=args.pop_refalt_af_col,
                               chrom_col="CHROM", pos_col="POS",
                               ref_col="REF", alt_col="ALT")
        # list of chromosomes to search for
        chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                  "11","12","13","14","15","16","17","18","19","20",
                  "21","22","X","Y"]
        # for each chrom ..
        for chr_i in chroms:
            # define chrom-level filepath
            pop_refalt_tsv_i = args.pop_refalt_prefix + ".chr" + chr_i + ".tsv.gz"

            # load chromosome-level tsv
            pop_refalt.load_refalt_tsv(pop_refalt_tsv_i)

    """
    init vars for tracking number of allele/strand flips
    """
    nvar=0
    nvar_match=0
    nvar_strandflip=0
    nvar_alleleflip=0
    nvar_strandalleleflip=0
    nvar_ambig_popfreq_match=0
    nvar_ambig_popfreq_alleleflip=0
    nvar_out=0

    """
    iterate through daner rows, and for each one, align A1 to reference allele,
    A2 to the alternate allele.
    """
    header_list = daner.header_list
    header_str = args.out_delim.join(header_list)
    print(header_str)
    while(1):
        daner.get_row(return_dict=True)
        if daner.row_str == "":
            break
        elif daner.row_str[0] == "#":
            continue

        # increment var count by 1
        nvar += 1

        # form instance of marker obj
        marker = Marker(daner.row_dict[args.varid_col],
                        chr=daner.row_dict[args.chr_col],
                        bp=daner.row_dict[args.pos_col],
                        a1=daner.row_dict[args.a1_col],
                        a2=daner.row_dict[args.a2_col],
                        eff=daner.row_dict[effsize_col],
                        eff_type=effsize_col)

        # add AFs to marker obj
        if args.ctrl_freq_col != None:
            try:
                marker.frq_u = float(daner.row_dict[args.ctrl_freq_col])
                marker.frq = float(daner.row_dict[args.ctrl_freq_col])
            except:
                sys.exit("\nERROR : non-numeric AFs detected. Ensure " + \
                         "all AFs are numeric before running.\n")
        if args.case_freq_col != None:
            try:
                marker.frq_a = float(daner.row_dict[args.case_freq_col])
            except:
                sys.exit("\nERROR : non-numeric AFs detected. Ensure " + \
                         "all AFs are numeric before running.\n")

        # only keep single nucleotide variants (ACGT). indels not supported for now.
        if marker.is_snv == False: continue

        # get coord info
        chrom = marker.chr
        pos = int(marker.bp)
        pos0 = pos - 1

        # adjust formatting of chrom based on whether or not ref contigs have 'chr'
        x_is_withchr = chrom.find("chr") == 0
        if is_withchr == True and x_is_withchr == False:
            chrom = 'chr' + chrom
        elif is_withchr == False and x_is_withchr == True:
            chrom = chrom.replace("chr","")

        # get ref allele at chr/pos, change to uppercase
        if is_2bit == True:
            ref = ref_fh[chrom].get_slice(pos0, pos)
        else:
            ref = ref_fh.fetch(chrom, pos0, pos)
        ref = ref.upper()

        # if a1a2 are ambiguous, cannot be certain that solution here is
        # accurate. Can only rescue if variant is in pop refalt tsv. flip
        # alleles if allele freqs in daner and in pop refalt tsv are flipped.
        
        if marker.ambiguous == False:
            # See if the reference allele is either a1 or a2, as listed.
            # If reference allele is not a1 or a2 as listed, try a1f or a2f.
            # If still no match, then call error.
            if ref == marker.a2:
                # no strand or allele flipping required
                nvar_match += 1
            elif ref == marker.a1:
                # allele flip is required
                marker.allele_flip()
                nvar_alleleflip += 1
            elif ref == marker.a2f:
                # strand flip is required
                marker.strand_flip()
                nvar_strandflip += 1
            elif ref == marker.a1f:
                # a strand flip followed by an allele flip required
                marker.strand_flip()
                marker.allele_flip()
                nvar_strandalleleflip += 1
            else:
                sys.exit("\nERROR : cannot map reference nucleotide at site to daner entry.\n")

        else:
            # skip if control freq column not defined
            if args.ctrl_freq_col == None:
                continue
            # skip if pop refalt not defined
            if pop_refalt == None:
                continue

            # define alt
            if ref == marker.a2:
                alt = marker.a1
            elif ref == marker.a1:
                alt = marker.a2
            else:
                sys.exit("\nERROR : cannot map reference nucleotide at site to daner entry.\n")

            # get refalt dir from pop refalt data 
            refalt_dir = pop_refalt.refalt_dir(marker.chr, 
                                               marker.bp, 
                                               ref, alt)
            # get a1 control freq in daner file
            a1_freq = float(daner.row_dict[args.ctrl_freq_col])

            # if ambiguous AF exclude range defined, skip if variant falls
            # within defined AF range
            if ambiguous_af_range_exclude != None:
                if (ambiguous_af_range_exclude[0] <= a1_freq) and (a1_freq <= ambiguous_af_range_exclude[1]):
                    continue

            # get a1 a2 pop direction
            a1a2_dir = a1_a2_direction(marker.a2, 
                                       marker.a1, 
                                       a1_freq)
            
            # only proceed if directions aren't set as 'none'
            if a1a2_dir == None or refalt_dir == None:
                continue

            # if refalt allelic balance doesn't match a1a2, then perform allele flip
            if a1a2_dir != refalt_dir:
                marker.allele_flip()
                nvar_ambig_popfreq_alleleflip += 1
            else:
                nvar_ambig_popfreq_match += 1

        # after strand/allele flips, ref and alt should now be a2 and a1 respectively

        """
        apply aligned A1/A2 and effect stats to row
        """
        daner.row_dict[args.a1_col] = marker.a1
        daner.row_dict[args.a2_col] = marker.a2
        daner.row_dict[effsize_col] = round(marker.eff, args.round)
        try:
            daner.row_dict[args.case_freq_col] = round(marker.frq_a, args.round)
        except:
            pass
        try:
            daner.row_dict[args.ctrl_freq_col] = round(marker.frq_u, args.round)
        except:
            pass

        """
        reconstruct row list
        """
        row_list = []
        for col in header_list:
            row_list.append(str(daner.row_dict[col]))
        row_str = args.out_delim.join(row_list)
        print(row_str)

        # inrement total output variant count by 1
        nvar_out += 1

    daner.close_fh()

    """
    if defined by user, write nvar stats to log file
    """
    if args.log_file != None:
        log_fh = open(args.log_file, "w")
        log_fh.write("nvar_input\t" + str(nvar) + "\n")
        log_fh.write("nvar_match\t" + str(nvar_match) + "\n") 
        log_fh.write("nvar_strandflip\t" + str(nvar_strandflip) + "\n")
        log_fh.write("nvar_alleleflip\t" + str(nvar_alleleflip) + "\n")
        log_fh.write("nvar_strandalleleflip\t" + str(nvar_strandalleleflip) + "\n")
        log_fh.write("nvar_ambig_popfreq_match\t" + \
                     str(nvar_ambig_popfreq_match) + "\n")
        log_fh.write("nvar_ambig_popfreq_alleleflip\t" + \
                     str(nvar_ambig_popfreq_alleleflip) + "\n")
        log_fh.write("nvar_output\t" + str(nvar_out) + "\n")
        log_fh.close()

    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--in-delim", type=str,
                      action="store", default=None, 
                      help="Delimiter for input daner.")
    args.add_argument("--out-delim", type=str,
                      action="store", default="\t")
    args.add_argument("--chr-col", type=str,
                      action="store", default="CHR")
    args.add_argument("--pos-col", type=str,
                      action="store", default="BP")
    args.add_argument("--varid-col", type=str,
                      action="store", default="SNP")
    args.add_argument("--a1-col", type=str,
                      action="store", default="A1")
    args.add_argument("--a2-col", type=str,
                      action="store", default="A2")
    args.add_argument("--effsize-col", type=str,
                      action="store", default="BETA")
    args.add_argument("--case-freq-col", type=str,
                      action="store", default=None,
                      help="column name for case a1 freq in " + \
                           "input daner file.")
    args.add_argument("--ctrl-freq-col", "--freq-col", type=str,
                      action="store", default=None,
                      help="column name for control or overall a1 freq in " + \
                           "input daner file. Required for pop freq-based " + \
                           "rescue of ambiguous variants. If outcome is " + \
                           "continuous, place column for overall frequency here.")
    args.add_argument("--round", type=int, 
                      action="store", default=6,
                      help="number of places to round effsize values to, if flipped.")
    args.add_argument("--ambiguous-af-range-exclude", action="store", type=str, default=None,
                      help="For ambiguous (A/T, G/C) SNVs, allele frequency " + \
                           "range in daner file to subset on, " + \
                           "in form of LOWERLIMIT-UPPERLIMIT. " + \
                           "If range is defined, an ambiguous SNV must fall " + \
                           "outside of range in daner file to be considered " + \
                           "for inclusion.")
    args.add_argument("--pop-refalt-prefix", type=str,
                      action="store",
                      default=None,
                      help="prefix name for " + \
                           "tab-delimited files with 1000 genomes " + \
                           "phase 3 simplified allele frequencies. For a " + \
                           "variant in a given population, '2' means " + \
                           "reference allele is more common, '1' means " + \
                           "alternate allele is more common, '0' means " + \
                           "ref and alt allele are equally common [default %(default)s]")
    args.add_argument("--pop-refalt-af-col", type=str,
                      action="store", default="EUR_AF",
                      help="column from pop refalt tsv to use for allele " +\
                           "frequency based alignment")
    args.add_argument("--log-file", type=str, action="store", default=None,
                      help="log file to write allele/strand flip statistics to.")
    args.add_argument("ref_file", type=str,
                      help="reference genome file. The following filetypes are "+ \
                           "supported  : fasta, fasta.gz, fa, fa.gz, 2bit")
    args.add_argument("daner_file", type=str,
                      help="GWAS summary statistics file in the daner file format.")
    
    return args.parse_args()

if __name__ == "__main__":
    main()

