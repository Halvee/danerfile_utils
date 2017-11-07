import gzip
import operator

class Tbl(object):
    def __init__(self, filename, with_header=True, 
                 cols_recode_str=None):
        self.filename = filename
        self.fh = None
        self.header_list = None
        self.col_idx = {}
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
            self.fh = gzip.open(self.filename, "rb")
        else:
            self.fh = open(self.filename, "r")
        return self

    def header_to_idx(self, delim=None):
        header_str = self.fh.readline().rstrip()
        if delim != None:
            self.header_list = header_str.split(delim)
        else:
            self.header_list = header_str.split()
        for i in range(len(self.header_list)):
            if self.header_list[i] in self.cols_recode_dict:
                col_old = self.header_list[i]
                print(self.header_list[i])
                self.header_list[i] = self.cols_recode_dict[col_old]
                print(self.header_list[i])
            self.col_idx[ self.header_list[i] ] = i
        return self

    def get_row(self, delim=None, return_dict=False,
                return_list_too=False):
        row_str = self.fh.readline().rstrip()
        if delim != None:
            row_list = row_str.split(delim)
        else:
            row_list = row_str.split()
        if return_dict == False:
            return row_list
        else:
            row_dict = {}
            if len(self.col_idx) != len(row_list):
                return row_dict
            for col_name in self.col_idx:
                i = self.col_idx[col_name]
                row_dict[col_name] = row_list[i]
            if return_list_too == True:
                return row_dict, row_list
            else:
                return row_dict

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
            try:
                data = line.rstrip().split()
                [param, operand_str, value] = data[:3]
            except:
                raise Exception("Improper formatting in " + self.cnds_file + ": " + \
                                "line " + str(i) + " : " + line.rstrip())
            if operand_str in ["in","nin"]:
                operand = operand_str
                if os.path.isfile(value):
                    value = read_set_file(value)
                else:
                    value = value.split(",")
            else:
                print(operator.__dict__[operand_str])
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

    def test(self, val_dict):
        for cnd in self.cnds:
            [param, operand, value] = cnd
            val = val_dict[param]
            if operand == "in":
                if val not in value:
                    return False
            elif operand == "nin":
                if val in value:
                    return False
            else:
                if  val.replace(".","").isdigit():
                    if operand(float(val), float(value)) == False:
                        return False
                else:
                    if operand(val, value) == False:
                        return False
        return True

