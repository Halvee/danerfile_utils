import sys
import argparse
from lib.tbl import Tbl 
from lib.daner import Marker
try:
    from scipy.stats import binom_test
except:
    print("ERROR : required module 'scipy' not installed.")
    sys.exit(1)

"""
REQUIRED INPUT
 1. A list of clumped variants from your source case/control GWAS.
    Set up for clumped file produced by PLINK. Expected to have a 
    min(P) column, corresponding to the SNP in the clump with the lowest
    P-value.
 2. A PGC daner-formatted summary stats file for your target case/control
    GWAS. Standard file format expected.
       
WORKFLOW
 1. Read clumped source case/control GWAS file, retaining only the clumps
    where min(P) < args.max_p_clump. Store corresponding odds ratio or 
    regression beta coefficient.
 2. Read PGC daner-formatted summary stats file, retaining only the variants
    (chr-pos-ref-alt) that are min(P) from a clump from source file. Store
    corresponding odds ratio or regression beta coefficient.
 3. Get sum of vars that have same case/control direction of effect in 
    source and target files. Null hypothesis == 50% of vars. Binomial 
    test used to see if data supports null hypothesis.
"""
def main():
    args = parse_args()

    """ 
    define max p-value threshold from user input.
    """
    try:
        pval_thresholds = args.pval_thresholds.split(",")
        pval_thresholds = [float(pval) for pval in pval_thresholds]
    except:
        print("ERROR : input '" + args.pval_thresholds + \
              "' not in the form of comma delimited p-value threshold " + \
              "float values.")
        sys.exit(1)
    pval_thresholds.sort()
    pval_max = max(pval_thresholds)

    """
    create instance of Tbl class for source cc clumped file. Get p-value
    per clump, associated marker ID tagging min(P).
    """
    varids_old_keep = set()
    src_clump = Tbl(args.source_clumped_file, with_header=True)
    while(1):
        row = src_clump.get_row(delim=None, return_dict=True)
        if len(row) == 0: 
            break
        if float(row["P"]) < pval_max:
            varids_old_keep.add(row["SNP"])
        
    src_clump.close_fh()

    """ 
    create Tbl instance for source cc marker-level stats. read in data.
    create dictionary of case/control effect sizes.
    """
    src_daner = Tbl(args.source_daner_file, with_header=True)
    src_data = {}
    while(1):
        row = src_daner.get_row(delim=None, return_dict=True)
        if len(row) == 0:
            break
        varid_old = row["SNP"]
        if varid_old not in varids_old_keep:
            continue

        # NOTE : row chr/bp/a1/a2 should match clump info for marker
        src_data[varid_old] = Marker(row["SNP"],
                                     chr=row["CHR"], bp=int(row["BP"]), 
                                     a1=row["A1"], a2=row["A2"], 
                                     eff=float(row[args.eff_type]), 
                                     eff_type=args.eff_type,
                                     p=float(row["P"]))

    src_daner.close_fh()

    """
    from source data, create new dictionary where chr-bp-a1-a2 are keys
    """
    src_data_t = dict()
    for varid_old in src_data:
        var_data = src_data[varid_old]
        varid_new = "-".join([var_data.chr,
                              str(var_data.bp),
                              var_data.a1,
                              var_data.a2])
        src_data_t[varid_new] = var_data

    """
    create Tbl instance for target cc marker-level stats. read in data.
    """
    trg_daner = Tbl(args.target_daner_file, with_header=True)
    markers_cmp = []
    while(1):
        row = trg_daner.get_row(delim=None, return_dict=True)
        if len(row) == 0:
            break
        varid_a1a2 = "-".join([row["CHR"], row["BP"], row["A1"], row["A2"]])
        varid_a2a1 = "-".join([row["CHR"], row["BP"], row["A2"], row["A1"]])
        if varid_a1a2 not in src_data_t and varid_a2a1 not in src_data_t:
            continue
        
        marker_trg = Marker(row["SNP"],
                            chr=row["CHR"],
                            bp=int(row["BP"]),
                            a1=row["A1"],
                            a2=row["A2"],
                            p=float(row["P"]),
                            eff=float(row[args.eff_type]),
                            eff_type=args.eff_type)
        
        cmp_pass = False
        if marker_trg.name_ref in src_data_t:
            cmp_pass = True
        else:
            marker_trg.allele_flip()
            if marker_trg.name_ref in src_data_t:
                cmp_pass = True
        
        if cmp_pass == True:
            marker_src = src_data_t[marker_trg.name_ref]
            if marker_trg.eff_dir == "=" or marker_src.eff_dir == "=":
                continue
            markers_cmp.append([marker_src.name_ref, 
                                marker_src.eff_dir,
                                marker_trg.eff_dir])


    trg_daner.close_fh()

    """
    count the number of markers where the direction of effect is the same.
    """
    n = len(markers_cmp)
    n_samedir = 0
    for marker in markers_cmp:
        if marker[1] == marker[2]:
            n_samedir += 1

    """
    perform binomial test on obs vs exp frac samedir
    """
    binom_p = binom_test(n_samedir,n=n,p=0.5)
    print("n = "  + str(n))
    print("n(same direction) = " + str(n_samedir))
    print("binom(p) = " + str(binom_p))
    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--pval-thresholds", type=str,
                      action="store", default="0.0001", 
                      help="Comma-delimited list of p-val thresholds to eval.")
    args.add_argument("--eff-type", type=str, choices=("OR","BETA"),
                      action="store", default="OR",
                      help="Effect direction type to use, corresponds to col names.")
    args.add_argument("source_clumped_file", 
                      help="Source LD clumped file for case/control, made by PLINK, " + \
                           "with min(P) variant tagged.")
    args.add_argument("source_daner_file",
                      help="Source case/control daner file with odds ratio or beta value column included.")
    args.add_argument("target_daner_file",
                      help="Target case/control daner file with odds ratio or beta value column included.")
    
    return args.parse_args()

def eff_dir(eff, eff_type="OR", flipped=False):
    eq_val = 1
    if eff_type == "BETA":
        eq_val = 0
    if flipped == True:
        if eff < eq_val: return "-"
        elif eff > eq_val: return "+"
        else: return "="
    else:
        if eff < eq_val: return "+"
        elif eff > eq_val: return "-"
        else: return "="
 

if __name__ == "__main__":
    main()

