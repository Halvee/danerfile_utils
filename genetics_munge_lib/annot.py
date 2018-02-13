
import sys

IMPACT_RANKINGS = {"HIGH":3,
                   "MODERATE":2,
                   "LOW":1,
                   "MODIFIER":0,
                   "NA":-1}
IMPACT_RANKINGS_INV = {3:"HIGH",
                       2:"MODERATE",
                       1:"LOW",
                       0:"MODIFIER",
                       -1:"NA"}

class AnnotTxs(object):
    def __init__(self, annot_keys, annot_line, 
                 delim=",", subdelim="|", format="ANN"):
        self.annot_keys = annot_keys
        self.annot_line = annot_line
        self.delim = delim
        self.subdelim = subdelim
        self.format = format
        self.max_eff = None
        self.load_annot_line()
    
    def load_annot_line(self):
        self.annots = self.annot_line.split(self.delim)
        for i in range(len(self.annots)):
            tx = self.annots[i]
            self.annots[i] = AnnotTx(self.annot_keys, tx, delim=self.subdelim)
        return self

    def max_csq(self, impact_max=True, 
                sift_min=False, polyphen_max=False):
        global IMPACT_RANKINGS_INV
        annot_classifs = {}
        max_annot_ranking = -1
        max_gene_symbols_set = set()
        max_ppn = 0
        min_sift = 1
        for i in range(len(self.annots)):
            annot_impact = self.annots[i].__dict__["IMPACT"]
            if annot_impact == "": annot_impact = "NA"
            impact_ranking = IMPACT_RANKINGS[annot_impact]
            if impact_ranking not in annot_classifs: 
                annot_classifs[impact_ranking] = set()
            annot_classifs[impact_ranking].add(i)
            if sift_min == True:
                annot_sift = self.annots[i].__dict__["SIFT"]
                if annot_sift == "": annot_sift = "1"
                annot_sift = float(annot_sift)
                if annot_sift < min_sift: min_sift = annot_sift
            if polyphen_max == True:
                annot_polyphen = self.annots[i].__dict__["PolyPhen"]
                if annot_polyphen == "": annot_polyphen = "0"
                annot_polyphen = float(annot_polyphen)
                if annot_polyphen > max_ppn: max_ppn = annot_polyphen

        
        if len(annot_classifs) == 0:
            max_annot_ranking = -1
            max_gene_symbols="NA"
        else:
            max_annot_ranking = max(annot_classifs.keys())
            for i in annot_classifs[max_annot_ranking]:
                annot = self.annots[i]
                max_gene_symbols_set.add(annot.SYMBOL)
            if "" in max_gene_symbols_set: max_gene_symbols_set.remove("")
            max_gene_symbols = list(max_gene_symbols_set)
            max_gene_symbols.sort()
            max_gene_symbols = ",".join(max_gene_symbols)
            if max_gene_symbols == "": max_gene_symbols = "NA"

        max_annot_ranking_classif = IMPACT_RANKINGS_INV[max_annot_ranking]
        max_out = {"CSQ_MAX_IMPACT": max_annot_ranking_classif,
                   "CSQ_MAX_IMPACT_GENESYMBOL": max_gene_symbols,
                   "CSQ_MAX_PolyPhen": max_ppn,
                   "CSQ_MAX_SIFT": min_sift}

        return max_out

class AnnotTx(object):
    def __init__(self, keyvals, 
                  vals_str, delim="|"):
        for keyval in keyvals:
            self.__dict__[keyval] = None
        if delim != None:
            vals = vals_str.split(delim)
        else:
            vals = vals_str
        if len(keyvals) != len(vals):
            raise Exception("n(annot keys) != n(annot vals)")
        for i in range(len(vals)):
            key_i = keyvals[i]
            val_i = vals[i]
            self.__dict__[key_i] = val_i

def process_csq_header_desc(header_desc_str,                                    
                            replace_chars=["'",'"'," "],                        
                            main_split="Format: ",                              
                            subsplit="|"):                                      
    header_desc_list = header_desc_str.split(main_split)                        
    header_desc_str = header_desc_list[-1]                                      
    for replace_char in replace_chars:                                          
        header_desc_str = header_desc_str.replace(replace_char, "")             
    header_desc_list = header_desc_str.split(subsplit)                          
    return header_desc_list     
