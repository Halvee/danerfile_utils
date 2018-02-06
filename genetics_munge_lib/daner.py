import gzip
import operator

class Marker(object):
    def __init__(self, name, 
                 chr=None, bp=None, 
                 a1=None, a2=None,
                 eff=None, eff_type="OR",
                 p=None, se=None, info=None,
                 line=None):
        self.name = name
        self.chr = chr
        self.bp = int(bp)
        self.a1 = a1
        self.a2 = a2
        self.name_ref = "-".join([chr,str(bp),a1,a2]) 
        self.eff_type = eff_type
        self.eff = eff
        self.se = se
        self.p = p
        self.info = info
        self.line = line
        self.eff_dir = None
        self.get_eff_dir()

    def allele_flip(self):
        a1_old = self.a1
        a2_old = self.a2
        self.a1 = self.a2_old
        self.a2 = self.a1_old
        self.name_ref = "-".join([self.chr,self.bp,self.a1,self.a2])
        if self.eff_type == "BETA":
            self.eff = self.eff * -1
        else:
            self.eff = 1.0 / self.eff
        self.get_eff_dir()
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
