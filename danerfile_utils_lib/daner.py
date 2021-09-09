import sys
import gzip
import operator

## global vars
STRANDFLIP_NTS={"T":"A",
                "A":"T",
                "G":"C",
                "C":"G"}
AMBIGUOUS_A1A2=set(["AT","TA","GC","CG"])

class Marker(object):
    def __init__(self, name, 
                 chr=None, bp=None, 
                 a1=None, a2=None,
                 eff=None, eff_type="OR",
                 frq=None, frq_a=None, frq_u=None,
                 p=None, se=None, info=None,
                 line=None):
        self.name = name
        self.chr = chr
        self.bp = int(bp)
        self.a1 = a1
        self.a2 = a2
        self.a1a2 = self.a1 + self.a2
        self.name_ref = "-".join([chr,str(bp),a1,a2])
        self.frq = frq
        self.frq_a = frq_a
        self.frq_u = frq_u 
        self.eff_type = eff_type
        try:
            self.eff = float(eff)
        except:
            self.eff = None
        self.se = se
        self.p = p
        self.info = info
        self.line = line
        self.eff_dir = None
        self.get_eff_dir()

        # get globals
        global STRANDFLIP_NTS
        global AMBIGUOUS_A1A2

        # if variant is standard snv, label as such
        self.is_snv = False
        if self.a1 in STRANDFLIP_NTS and self.a2 in STRANDFLIP_NTS:
            self.is_snv = True

        # determine if variant is ambiguous (A/T, G/C)
        self.ambiguous = False
        if self.a1a2 in AMBIGUOUS_A1A2:
            self.ambiguous = True

        # get strandflip alleles. only proceed if a1 and a2 are ATGC
        self.a1f = None
        self.a2f = None
        self.a1fa2f = None
        if self.a1 in STRANDFLIP_NTS and self.a2 in STRANDFLIP_NTS:
            self.a1f = STRANDFLIP_NTS[self.a1]
            self.a2f = STRANDFLIP_NTS[self.a2]
            self.a1fa2f = self.a1f + self.a2f

    def allele_flip(self):
        a1_old = self.a1
        a2_old = self.a2
        self.a1 = a2_old
        self.a2 = a1_old
        self.name_ref = "-".join([self.chr,str(self.bp),self.a1,self.a2])
        
        # numeric effect must be provided, if not, then halt
        try:
            assert isinstance(self.eff, float)
        except:
            sys.exit("ERROR : non-numeric effect sizes detected, " + \
                     "remove from input file and then rerun for proper allele flipping.")

        # flip effect directions
        if self.eff_type == "BETA":
            self.eff = self.eff * -1
        else:
            self.eff = 1.0 / self.eff
        self.get_eff_dir()
        # flip frequencies, if defined
        if self.frq != None: self.frq = 1 - self.frq
        if self.frq_a != None: self.frq_a = 1 - self.frq_a
        if self.frq_u != None: self.frq_u = 1 - self.frq_u
        return self

    def strand_flip(self):
        # dict of flipped nts
        global STRANDFLIP_NTS
        # only proceed if all A1/A2 chars are A/T/G/C
        a1a2 = self.a1 + self.a2
        for i in range(len(a1a2)):
            if a1a2[i] not in STRANDFLIP_NTS:
                raise ValueError('non-ATGC character in A1/A2 (' + self.a1 + \
                                 '/' + self.a2 + ')')
        # flip a1
        for i in range(len(self.a1)):
            self.a1[i] = STRANDFLIP_NTS[self.a1[i]]
        # flip a2
        for i in range(len(self.a2)):
            self.a2[i] = STRANDFLIP_NTS[self.a2[i]]
        return self

    def get_eff_dir(self):
        if self.eff == None: return self
        eq_val = 1
        if self.eff_type == "BETA":
            eq_val = 0
        if self.eff < eq_val:
            self.eff_dir = "-"
        elif self.eff > eq_val:
            self.eff_dir = "+"
        else:
            self.eff_dir = "="
        return self

    def strand_align(self, Marker_aln):
        chrbp_x = (self.chr, self.bp)
        a1a2_x = set([self.a1, self.a2])
        chrbp_y = (Marker_aln.chr, Marker_aln.bp)
        a1a2_y = set([Marker_aln.a1, Marker_aln.a2])
        if chrbp_x!=chrbp_y or a1a2_x!=a1a2_y:
            raise Exception("incompatible markers")

        if self.chr != Marker_aln.chr: return self
        elif self.bp != Marker_aln.bp: return self
        if self.a1 == Marker_aln.a1 and self.a2 == Marker_aln.a2:
            return self
        elif self.a1 == Marker_aln.a2 and self.a2 == Marker_aln.a1:
            self.allele_flip()
