import io
import os
import sys
import gzip
import operator
import misc

class Tbl(object):
    def __init__(self, filename, 
                 delim = None,
                 with_header=True, 
                 cols_recode_str=None):
        self.filename = filename
        self.fh = None
        self.header_str = None
        self.header_list = None
        self.delim=delim
        self.col_idx = {}
        self.row_str = None
        self.row_list = []
        self.row_dict = {}
        self.open_fh()
        if with_header == True: 
            self.load_cols_recode(cols_recode_str)
            self.header_to_idx()

    def load_cols_recode(self, cols_recode_str):
        self.cols_recode_dict = {}
        if cols_recode_str == None: return self
        cols_recode_list = cols_recode_str.split(",")
        for col_oldnew_str in cols_recode_list:
            (col_old,col_new) = col_oldnew_str.split(":")[:2]
            self.cols_recode_dict[col_old] = col_new
        return self

    def open_fh(self):
        if self.filename == "stdin":
            self.fh = sys.stdin
        elif self.filename.find(".gz") != -1:
            fh = gzip.open(self.filename, "rb")
            self.fh = io.BufferedReader(fh)
        else:
            self.fh = open(self.filename, "r")
        return self

    def header_to_idx(self, delim=None):
        self.header_str = self.fh.readline().rstrip()
        self.header_list = self.header_str.split(self.delim)
        for i in range(len(self.header_list)):
            if self.header_list[i] in self.cols_recode_dict:
                col_old = self.header_list[i]
                self.header_list[i] = self.cols_recode_dict[col_old]
            self.col_idx[ self.header_list[i] ] = i
        return self

    def get_row(self, 
                return_dict=False):
        self.row_str = self.fh.readline().rstrip()
        if self.row_str == "":
            return self
        elif self.row_str[0] == "#": 
            return self
        
        self.row_list = self.row_str.split(self.delim)
        if return_dict == True:
            self.row_dict = {}
            for col_name in self.col_idx:
                i = self.col_idx[col_name]
                self.row_dict[col_name] = self.row_list[i]
        return self

    def close_fh(self):
        self.fh.close()
        return self

class Cnds(object):
    def __init__(self, cnds_file):
        self.cnds_file = cnds_file
        self.cnds = []
        self.read_cnds_file()

    def read_cnds_file(self, print_to_stdout=False):
        if self.cnds_file == None: return self
        fh = open(self.cnds_file, "r")
        i = 0
        for line in fh:
            i += 1
            if line[0] == "#": continue
            try:
                data = line.rstrip().split()
                [param, operand_str, value] = data[:3]
                param = param.replace('"','')
                value = value.replace('"','')
            except:
                raise Exception("Improper formatting in " + self.cnds_file + ": " + \
                                "line " + str(i) + " : " + line.rstrip())
            if operand_str in ["in","nin"]:
                operand = operand_str
                if os.path.isfile(value):
                    value = read_set_file(value)
                else:
                    value = set(value.split(","))
            elif operand_str in ["grep","grepv"]:
                operand = operand_str
                value = value
            else:
                try:
                    operand = operator.__dict__[operand_str]
                except:
                    raise Exception("Improper formatting in " + self.cnds_file + ": " + \
                                    "operand " + operand_str + " not defined.")
            cnd = [param, operand, value]
            self.cnds.append(cnd)
            if print_to_stdout == True:
                out_str = " ".join(["CONDITION",":",param, operand_str, value])
                print(out_str)
                
        fh.close()
        return self

    def test(self, val_dict, forgive_missing_cnd=False):
        for cnd in self.cnds:
            [param, operand, value] = cnd
            if param not in val_dict and forgive_missing_cnd==True: continue
            try:
                val = val_dict[param]
            except:
                if forgive_missing_cnd == True:
                    continue
                else:
                    print("ERROR : cnds file param '"+param+"' missing from table.")
                    sys.exit(1)

            if operand == "in":
                if val not in value:
                    return False
            elif operand == "nin":
                if val in value:
                    return False
            elif operand == "eq":
                if operand(val, value) == False:
                    return False
            elif operand == "grep":
                if val.find(value) == -1:
                    return False
            elif operand == "grepv":
                if val.find(value) != -1:
                    return False
 
            else:
                if val.replace(".","").isdigit() and value.replace(".","").isdigit():
                    if operand(float(val), float(value)) == False:
                        return False
                else:
                    if operand(val, value) == False:
                        return False
        return True

def read_set_file(set_filename):
    set_out = set()
    fh = misc.open_file(set_filename)
    for line in fh:
        set_out.add(line.rstrip())
    fh.close()
    return set_out

