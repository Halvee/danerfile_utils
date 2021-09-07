

from .tbl import Tbl


class PopRefAlt(object):
    def __init__(self,
                 delim="\t",
                 pop_refalt_af_col="EUR_AF",
                 chrom_col="CHROM", pos_col="POS",
                 ref_col="REF", alt_col="ALT"):

        # init af tree
        self.pop_refalt = dict()

        # store delimiter for reading input tables
        self.delim  = delim

        # store column name for chrom/pos/ref/alt,
        # allele frequency
        self.chrom_col=chrom_col
        self.pos_col = pos_col
        self.ref_col = ref_col
        self.alt_col = alt_col
        self.af_col = pop_refalt_af_col

    def load_refalt_tsv(self,
                        pop_refalt_tsv):

        # init instance of Tbl
        pop_refalt_tbl = Tbl(pop_refalt_tsv,
                             delim=self.delim,
                             with_header=True)

        # main loop for reading pop refalt tsv
        while(1):
            pop_refalt_tbl.get_row(return_dict=True)
            if pop_refalt_tbl.row_str == "":
                 break
            elif pop_refalt_tbl.row_str[0] == "#":
                continue
            chrom = pop_refalt_tbl.row_dict[self.chrom_col]
            pos = int(pop_refalt_tbl.row_dict[self.pos_col])
            ref = pop_refalt_tbl.row_dict[self.ref_col]
            alt = pop_refalt_tbl.row_dict[self.alt_col]

            # strip 'chr' from chrom
            chrom = chrom.replace("chr","")

            # store af data to tree struct,
            # as follows :
            # 1 : alternate allele more common in the population
            # 2 : reference allele more common in the population
            # 0 : reference allele freq == alternate allele freq
            if chrom not in self.pop_refalt:
                self.pop_refalt[chrom] = dict()
            if pos not in self.pop_refalt[chrom]:
                self.pop_refalt[chrom][pos] = dict()
            if ref not in self.pop_refalt[chrom][pos]:
                self.pop_refalt[chrom][pos][ref] = dict()
            if alt not in self.pop_refalt[chrom][pos][ref]:
                af = float(pop_refalt_tbl.row_dict[self.af_col])
                if af < 0.5:
                    val = 2
                elif af > 0.5:
                    val = 1
                else:
                    val = 0
                self.pop_refalt[chrom][pos][ref][alt] = val

        return self

    def refalt_dir(self, chrom, pos, ref, alt):
        chrom = chrom.replace('chr','')
        if chrom not in self.pop_refalt: 
            return None
        if pos not in self.pop_refalt[chrom]: 
            return None
        if ref not in self.pop_refalt[chrom][pos]:
            raise ValueError("At site "+str(chrom)+":"+str(pos)+", " + \
                             "reference allele " + ref + " not found.")
        if alt not in self.pop_refalt[chrom][pos][ref]: 
            return None
        
        # retrieve refalt direction :
        # 2 means reference allele is more common
        # 1 means alternate allele is more common
        # 0 means ref and alt allele are equally common
        refalt_dir = self.pop_refalt[chrom][pos][ref][alt]
        return refalt_dir

def a1_a2_direction(a1, a2, a1_freq, 
                    undef_lowerbound=0.4, 
                    undef_upperbound=0.6):
    if a1_freq < undef_lowerbound:
        return 2
    elif a1_freq > undef_upperbound:
        return 1
    else:
        return 0
