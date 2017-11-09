
import sys

IMPACT_RANKINGS = {"HIGH":3,
                   "MODERATE":2,
                   "LOW":1,
                   "MODIFIER":0,
                   "NA":-1}

class SnpeffAnnotTxs(object):
    def __init__(self, snpeff_annot_line, delim=",", subdelim="|", format="ANN"):
        self.delim = delim
        self.subdelim = subdelim
        self.format = format
        self.max_eff = None
        self.load_snpeff_annot_line(snpeff_annot_line)
    
    def load_snpeff_annot_line(self, snpeff_annot_line):
        if self.format == "ANN":
            self.snpeff_annots = snpeff_annot_line.split(self.delim)
       	    for i in range(len(self.snpeff_annots)):
                snpeff_tx = self.snpeff_annots[i].split(self.subdelim)
                self.snpeff_annots[i] = SnpeffAnnotTx(snpeff_tx, 
                                                      format=self.format)
        return self

    def get_max_impact(self, protein_coding_only=False):
        global IMPACT_RANKINGS
        max_impact = "NA"
        for i in range(len(self.snpeff_annots)):
            if protein_coding_only == True:
                if self.snpeff_annots[i].tx_biotype!="protein_coding":
                    continue
            if self.snpeff_annots[i].impact == None: continue
            rank_i = IMPACT_RANKINGS[self.snpeff_annots[i].impact]
            rank_max = IMPACT_RANKINGS[max_impact]
            if rank_i > rank_max:
                max_impact = self.snpeff_annots[i].impact
                self.max_eff = self.snpeff_annots[i]
        return self

class SnpeffAnnotTx(object):
    def __init__(self, tx_l, format="ANN"):
        if format=="ANN":
            try:
                [self.alt,
                self.eff,
                self.impact,
                self.gene_name,
                self.gene_id,
                self.feature_type,
                self.feature_id,
                self.tx_biotype,
                self.exon_total,
                self.hgvs_c,
                self.hgvs_p,
                self.cdna_pos,
                self.cds_pos,
                self.prot_pos] = tx_l[:14]
            except:
                self.alt=None
                self.eff=None
                self.impact=None
                self.gene_name=None
                self.gene_id=None
                self.feature_id=None
                self.tx_biotype=None
                self.exon_total=None
                self.hgvs_c=None
                self.hgvs_p=None
                self.cdna_pos=None
                self.cds_pos=None
                self.prot_pos=None

