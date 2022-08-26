import os
import numpy as np
import pandas as pd



class IndTable:
    """
    Class to manipulate the indices table.

    Args:
        table_csv (str, None): Path to the csv table with indices data. If ``None`` use built-in table (default).

    *Class attributes:*

    Attributes:
        table (pd.DataFrame): Table with indices.
        indices (list): List of available indices.
        params (list): List of required index line parameters.
    """

    def __init__(self, table_csv=None, verb=False):

        if not table_csv:
            table_csv = os.path.join(os.path.dirname(__file__), "actin_table.csv")

        self.table = pd.read_csv(table_csv, index_col=False)
        self.indices = list(np.unique(self.table.ind_id))
        self.params = list(self.table.keys())


    def add_line(self, ind_id, ind_var, ln_id, ln_c, ln_ctr, ln_win, bandtype):
        """Add a line to the indices table.

        Args:
            ind_id (str): Index identification for the line.
            ind_var (str): A string identifying the line as activity line (if ``L1``, ``L2``, ...) or reference band (if ``R1``, ``R2``, ...).
            ln_id (str): Line identification
            ln_c (float): Constant value to be multiplied to the line.
            ln_ctr (float): Line centre [angstrom]
            ln_win (float): Line bandwidth [angstrom]
            bandtype (str): Bandpass type. If ``sq`` is square filter, if ``tri`` is triangular filter where ``ln_win`` will be the triangle FWHM.
        """
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
