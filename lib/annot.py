
IMPACT_RANKINGS = {"HIGH":3,
                   "MODERATE":2,
                   "LOW":1,
                   "MODIFIER":0}

class SnpeffAnnotTxs(object):
    def __init__(self, snpeff_annot_line, delim=",", subdelim="|", format="ANN"):
        self.delim = delim
        self.subdelim = subdelim
        self.format = format
        self.load_snpeff_annot_line(snpeff_annot_line)
        self.snpeff_annots = snpeff_annot_line.split(delim)   
        self.max_eff = None
        for i in range(len(self.snpeff_annots)):
            self.snpeff_annots[i] = self.snpeff_annots[i].split(subdelim)
    
    def load_snpeff_annot_line(self, snpeff_annot_line):
        if self.format == "ANN":
            self.snpeff_annots = snpeff_annot_line.split(self.delim)
       	    for i in range(len(self.snpeff_annots)):
                snpeff_tx = self.snpeff_annots[i].split(self.subdelim)
                self.snpeff_annots[i] = SnpeffAnnotTx(snpeff_tx, 
                                                      format=self.format)
        return self

    def get_max_impact(self):
        global IMPACT_RANKINGS
        max_impact = "MODIFIER"
        for i in range(len(self.snpeff_annots)):
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
                self.feature_id,
                self.tx_biotype,
                self.exon_total,
                self.hgvs_c,
                self.hgvs_p,
                self.cdna_pos,
                self.cds_pos,
                self.prot_pos,
                self.feature_dist,
                self.err] = tx_l
            except:
                print("ERROR : improperly formatted " + \
                      "snpeff ann entry for tx")
                print(tx_l)
                sys.exit(1)

