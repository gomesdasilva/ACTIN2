import os
import numpy as np
import pandas as pd



class IndTable:
    """
    Class to manipulate the indices table.

    Parameters:
    -----------
        table_csv : str
            Path to the csv table with indices data. If 'None' use built-in table.
    """

    def __init__(self, table_csv=None, verb=False):
        if verb:
            print("Running IndTable")

        if not table_csv:
            table_csv = os.path.join(os.path.dirname(__file__), "actin_table.csv")

        self.table = pd.read_csv(table_csv, index_col=False)
        
        self.indices = np.unique(self.table.ind_id)
        self.params = list(self.table.keys())

    def show_table(self):
        print(self.table)

    def show_indices(self):
        print(self.indices)

    def show_params(self):
        print(self.params)

    def add_line(self, ind_id, ind_var, ln_id, ln_c, ln_ctr, ln_win, bandtype):
        # parameters should be str or float
        new_line = dict(
            ind_id = ind_id,
            ind_var = ind_var,
            ln_id = ln_id,
            ln_c = ln_c,
            ln_ctr = ln_ctr,
            ln_win = ln_win,
            bandtype = bandtype
        )
        self.table = self.table.append(new_line, ignore_index=True)

    def del_line(self, ln_id):
        self.table = self.table[self.table["ln_id"] != ln_id]
        self.table.reset_index(inplace=True)

    def del_index(self, ind_id):
        self.table = self.table[self.table["ind_id"] != ind_id]
        self.table.reset_index(inplace=True)

    def save_table(self, filename):
        self.table.to_csv(filename, index=False)

    def get_index(self, index):
        return self.table[self.table.ind_id == index]
