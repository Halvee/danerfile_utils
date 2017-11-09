
import argparse
from lib import vcf,tbl,misc

def main():
    args = parse_args()
    
    printed_header = False

    cnds = tbl.Cnds(args.cnds_file)

    fam_trios = misc.load_fam_trios(args.fam_trios_file)

    for vcf_file in args.vcf_files:
        
        vcf_fh = misc.open_file(vcf_file)
        
        vcf_r = vcf.VcfReader(vcf_fh, cols_replace_format=args.cols_replace_format,
                              cols_replace_info=args.cols_replace_info)

        while(1):
            dat = vcf_r.next_line()
            if dat[1] == 0:
                pass
            elif dat[1] == 1:
                vcf_header = ["VARIANT_ID", "SAMPLE_ID"] + vcf_r.vcf_header[:7] + \
                              vcf_r.metainfo_lists["INFO"] + \
                              vcf_r.metainfo_lists["FORMAT"]
                full_header = vcf_header + ["TRANS_CLASSIF"]
                print(args.delim.join(full_header))
            elif dat[1] == 2:
                dat_rows = vcf_r.vcf_entry.get_sample_rows(vcf_r.metainfo_lists, vcf_r.sample_list)
                dat_rows_dict = {}
                for dat_row in dat_rows:
                    dat_dict = misc.keyval_list_pair_to_dict(vcf_header,dat_row)
                    dat_rows_dict[dat_dict["SAMPLE_ID"]] = dat_dict
                for proband in fam_trios:
                    proband_call = dat_rows_dict[proband]
                    father_call = dat_rows_dict[fam_trios[proband][0]]
                    mother_call = dat_rows_dict[fam_trios[proband][1]]
                    if cnds.test(proband_call) == False:
                        continue

                    # skip trios with no call evidence for variant
                    gt_list = [proband_call["GT"], father_call["GT"], mother_call["GT"]]
                    for i in range(len(gt_list)): 
                        gt_list[i] = gt_list[i].replace("/","|")
                        gt_list[i] = gt_list[i].replace(".","0")
                    gt_counts = get_gt_count(gt_list)
                    if sum(gt_counts) == 0 or sum(gt_counts) > 2: 
                        continue
                    if sum(gt_counts[1:]) > 1:
                        continue
                    
                    # skip var classif if min(dp) in trio is less than allowed
                    dp_trio = [proband_call["DP"], father_call["DP"], mother_call["DP"]]
                    for i in range(len(dp_trio)):
                        if dp_trio[i] == ".": dp_trio[i] = 0
                        else: dp_trio[i] = int(dp_trio[i])
                    if min(dp_trio) < args.min_dp_trio:
                        continue

                    # get ADs per sample
                    ad_trio = [proband_call["AD"], father_call["AD"], mother_call["AD"]]
                    multialleleic = False
                    for i in range(len(ad_trio)):
                        ad_trio[i] = ad_trio[i].split(",")
                        if len(ad_trio[i]) > 2: multialleleic = True
                        else:
                            ad_trio[i] = int(ad_trio[i][1])
                    if multialleleic == True: continue
                    
                    # get percent of reads carrying alt allele, get rid of trios where sum eq 0
                    alt_af_trio = []
                    for i in range(len(ad_trio)):
                        alt_af_trio.append( float(ad_trio[i]) / float(dp_trio[i]) )
                    if sum(alt_af_trio) == 0: continue

                    # get PL values per trio member
                    pl = [proband_call["PL"], father_call["PL"], mother_call["PL"]]
                    for i in range(len(pl)):
                        pl[i] = pl[i].split(",")
                        pl[i] = [int(j) for j in pl[i]]

                    # if father percent between 0.25 and 0.75..
                    trans_classif = "UNDEF"
                    perc_alt_range = [float(x) for x in args.het_perc_alt_range.split('-')]
                    
                    proband_dnm = pl[0][0] > args.pl_thresh and pl[0][1] < args.pl_thresh and pl[0][2] > args.pl_thresh
                    father_dnm = pl[1][0] < args.pl_thresh and pl[1][1] > args.pl_thresh and pl[1][2] > args.pl_thresh
                    mother_dnm = pl[2][0] < args.pl_thresh and pl[2][1] > args.pl_thresh and pl[2][2] > args.pl_thresh
                    
                    if proband_dnm and father_dnm and mother_dnm and misc.between(alt_af_trio[0], perc_alt_range):
                        if alt_af_trio[1] < args.dnm_perc_alt_min_parent and alt_af_trio[2] < args.dnm_perc_alt_min_parent:
                            trans_classif = "DENOVO"
                        elif alt_af_trio[2] == 0:
                            trans_classif = "MOSAIC_TRANS_FATHER"
                        elif alt_af_trio[1] == 0:
                            trans_classif = "MOSAIC_TRANS_MOTHER"
                        else:
                            trans_classif = "UNDEF"
                    elif misc.between(alt_af_trio[1], perc_alt_range):
                        if misc.between(alt_af_trio[0], perc_alt_range):
                            trans_classif = "HET_TRANS_FATHER"
                        else:
                            trans_classif = "HET_NTRANS_FATHER"
                    elif misc.between(alt_af_trio[2], perc_alt_range):
                        if misc.between(alt_af_trio[0], perc_alt_range):
                            trans_classif = "HET_TRANS_MOTHER"
                        else:
                            trans_classif = "HET_NTRANS_MOTHER"
                    else: 
                        trans_classif = "UNDEF"


                    #if trans_classif == "HIGH_CONF_MOSAIC":
                    #    print(dp_trio, ad_trio, alt_af_trio, proband_call["VARIANT_ID"])
                    #    print(proband_call["PL"], father_call["PL"], mother_call["PL"])
                    #    continue
                    if trans_classif == "UNDEF": continue
                    proband_call["TRANS_CLASSIF"] = trans_classif
                    out = []
                    for colname in full_header:
                        out.append(str(proband_call[colname]))
                    print("\t".join(out))
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

def get_gt_count(gts, delim="|"):
    if isinstance(gts, list):
        for i in range(len(gts)):
            gts[i] = sum([int(x) for x in gts[i].split(delim)])
    else:
        gts = sum([int(x) for x in gts.split(delim)])
    return gts

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
                      help="comma delimited set of old:new column names for FORMAT values.")
    args.add_argument("--min-dp-trio", type=int, action="store",
                      default=10, help="minimum allowed coverage in each trio member for a variant call.")
    args.add_argument("--cnds-file", type=str, action="store",
                      default=None, help="file with conditions to subset on.")
    args.add_argument('--het-perc-alt-range', type=str, action="store",
                      default="0.3-0.7", help="range for perc alt in reads to be het.")
    args.add_argument('--dnm-perc-alt-min-parent', type=float, action="store",
                      default=0.05, help="min perc alt allele in reads for parents, for DNMs.")
    args.add_argument('--pl-thresh', type=int, action="store", default=20,
                      help="thresholding to apply to PL.")
    args.add_argument("fam_trios_file", action="store", type=str,
                      help="name of fam file with trio relatedness info.")
    args.add_argument("vcf_files", nargs="+",
                      help="input VCF files.")
    
    return args.parse_args()

if __name__ == "__main__":
    main()
