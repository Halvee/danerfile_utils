
class VcfSite(object):
    def __init__(self, data_list):
        self.chrom = data_list[0]
        self.pos = int(data_list[1])
        self.id = data_list[2]
        self.ref = data_list[3]
        self.alt = data_list[4]

class VcfVariant(VcfSite):
    def __init__(self, data_list, cols_replace_info=None, **kwargs):
        VcfSite.__init__(self, data_list)
        try:
            self.qual = float(data_list[5])
        except:
            self.qual = 0
        self.filter = data_list[6]
        self.info = data_list[7]
        self.cols_replace_info = cols_replace_info
        self.info_dict = load_key_val(self.info,sep=";",subsep="=",
                                      replacements=load_key_val(self.cols_replace_info,
                                                                sep=",",
                                                                subsep=":"))
        self.info_list = list(self.info_dict.keys())
        self.info_list.sort()

class VcfGts(VcfVariant):
    def __init__(self, data_list, sample_list,
                 gt_param_delim=":", cols_replace_format=None,
                 **kwargs):
        VcfVariant.__init__(self, data_list, **kwargs)
        self.gt_param_delim=gt_param_delim
        self.format = data_list[8].split(self.gt_param_delim)
        self.cols_replace_format=cols_replace_format
        self.load_format_replacements()
        self.sample_gts = data_list[9:]
        self.gts = {}
        self.load_sample_gts(sample_list)

    def load_format_replacements(self):
        replacements=load_key_val(self.cols_replace_format,
                                  sep=",", subsep=":")
        if replacements == None:
            return self
        for i in range(len(self.format)):
            if self.format[i] in replacements:
                self.format[i] = replacements[self.format[i]]
        return self

    def load_sample_gts(self, sample_list):
        for i in range(len(sample_list)):
            sample_i = sample_list[i]
            self.gts[sample_i] = {}
            sample_i_vals = self.sample_gts[i].split(self.gt_param_delim)
            for j in range(len(self.format)):
                format_j = self.format[j]
                self.gts[sample_i][format_j] = sample_i_vals[j]
        return self

    def get_sample_rows(self, metainfo_lists, sample_list, varid_delim="-"):
        sample_rows = []
        var_row = [self.chrom, self.pos, 
                   self.id, self.ref, self.alt,
                   self.qual, self.filter]
        varid = varid_delim.join([self.chrom, str(self.pos),
                                  self.ref, self.alt])
        for col_name in metainfo_lists["INFO"]:
            if col_name not in self.info_list:
                var_row.append("0")
            else:
                var_row.append(self.info_dict[col_name])
        sample_rows = []
        for sampleid in sample_list:
            sample_row = [varid, sampleid] + var_row
            for col_name in metainfo_lists["FORMAT"]:
                if col_name not in self.gts[sampleid]:
                    sample_row.append("0")
                else:
                    sample_row.append(self.gts[sampleid][col_name])
            sample_rows.append(sample_row)

        return sample_rows


class VcfReader(object):
    def __init__(self, vcf_fh, delim="\t", 
                 cols_replace_info=None, 
                 cols_replace_format=None):
        self.vcf_fh = vcf_fh
        self.delim = delim
        self.metainfo = {}
        self.vcf_header = []
        self.cols_replace_info=cols_replace_info
        self.cols_replace_format=cols_replace_format
        self.cols_replace={"INFO":load_key_val(cols_replace_info,sep=",",subsep=":"),
                           "FORMAT":load_key_val(cols_replace_format,sep=",",subsep=":")}
        self.metainfo_lists = {"INFO":[],"FORMAT":[]}
        self.sample_list = []
        self.vcf_entry = None
        self.linenum = 0
        pass

    def load_metadata(self, line):
        info_start = line.find("<")
        info_end = line.find(">")
        metainfo_classif_end=line.find("=")
        metainfo_classif = line[:metainfo_classif_end]
        if metainfo_classif not in ("INFO","FORMAT"): 
            return self
        info = line[(info_start+1):info_end]
        keyval = load_key_val(info, sep=",", subsep="=")  
        if "ID" in keyval:
            if metainfo_classif not in self.metainfo:
                self.metainfo[metainfo_classif] = []
            self.metainfo[metainfo_classif].append(keyval)
            if self.cols_replace[metainfo_classif] != None:
                if keyval["ID"] in self.cols_replace[metainfo_classif]:
                    keyval["ID"] = self.cols_replace[metainfo_classif][keyval["ID"]]
            self.metainfo_lists[metainfo_classif].append(keyval["ID"])
        return self

    def load_header(self, line):
        self.vcf_header.extend(line.rstrip().split(self.delim))
        if len(self.vcf_header) > 9:
            self.sample_list = self.vcf_header[9:]
        return self

    def next_line(self):
        line = self.vcf_fh.readline()
        line = line.rstrip()
        if line == "":
            return self,3
        elif line[:2] == "##":
            self.load_metadata(line[2:])
            return self,0
        elif line[:1] == "#":
            self.load_header(line[1:])
            return self,1
        else:
            data = line.rstrip().split(self.delim)
            if len(data) > 0:
                self.vcf_entry = VcfGts(data, self.sample_list,
                                        cols_replace_info=self.cols_replace_info,
                                        cols_replace_format=self.cols_replace_format)
            elif len(data) > 5:
                self.vcf_entry = VcfVariant(data)
            else:
                self.vcf_entry = VcfSite(data)
            return self,2

def load_key_val(keyval_str, sep=";", subsep="=", replacements=None):
    if keyval_str == None: return None
    keyval = {}
    keyval_list = keyval_str.split(sep)
    for keyval_i in keyval_list:
        try:
            (key,val)=keyval_i.split(subsep)
            if replacements != None:
                key_new = key
                if key in replacements:
                    key_new = replacements[key]
            else:
                key_new = key
            keyval[key_new] = val
        except:
            keyval[keyval_i] = 1
    return keyval
