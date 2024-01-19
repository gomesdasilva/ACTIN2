import sys, os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import tqdm
import time

from .spectrographs._spec_tools import printif
from .calcindices import CalcIndices
from .readspec import ReadSpec
from .indtable import IndTable


class ACTIN:
    r"""The ACTIN class. Reads fits files and calculates activity indices.

    Any index is calculated as

    .. Math::
        I = \frac{\sum C_i F_i}{\sum K_j R_j}

    where :math:`F_i` and :math:`R_j` are the flux inside the activity sensitive lines and reference bands, respectively. :math:`C_i` and :math:`K_j` are constants to be multiplyed to the :math:`F_i` and :math:`R_j` fluxes.

    *Class attributes:*

    Attributes:
        ReadSpec (actin2.ReadSpec): Object that reads spectrum and   
            headers. See :py:class:`actin2.readspec.ReadSpec`
        IndTable (actin2.IndTable): Object containing the 
            indices parameters table. See :py:class:`actin2.indtable.IndTable`
        CalcIndices (actin2.CalcIndices): Object that calculates the 
            indices. See :py:class:`actin2.calcindices.CalcIndices`
    """

    def __init__(self):
        version_file = os.path.join(os.path.dirname(__file__), "VERSION")

        try:
            with open(version_file, 'r') as file:
                version = file.read()
        except FileNotFoundError:
            version = "unknown"

        self.version = version
        self.ReadSpec = ReadSpec
        self.IndTable = IndTable
        self.CalcIndices = CalcIndices


    def run(self, files, indices, obj_in=None, save_data=None, progress=True, verb=False, table_df=None, spec_kw=dict(), calcind_kw=dict()):
        """Run ACTIN for a list of fits files and indices.

        Args:
            files (list, str) : List of fits files paths, or string with one fits file path to be read.
            indices (list, str) : List of indices identification (as in 
                the indices table) to be calculated. If empty returns the fits headers only.
            rv_in (float, None) : If float, use this RV value to 
                shift spectrum to stellar rest frame, if None use RV value from fits file.
            spec_class_in (None, str) : Name of the spectrograph to be read. Usefull if new spectrographs are added to the ACTIN spectrograph directory.
            obj_in (str, None) : If not None, the target name from 
                the fits file is overriden by this input.
            save_data (str, None): path with file_name to save output in csv format.
            progress (bool) : If True show progress bar.
            table_df (None, pd.DataFrame) : Table DataFrame with the indices parameters.
                If None use the default table.
            verb (bool) : Activate verbose option.
        
        Returns:
            pandas.DataFrame : Output pandas table.
        """
        # If files is given as a string enforce type as list:
        if not isinstance(files, (list, np.ndarray, tuple)):
            files = [files]

        for file in files:
            if not os.path.isfile(file):
                raise FileNotFoundError(f"{__class__.__name__} ERROR: File '{file}' not found")

        if len(files) == 0:
            raise RuntimeError(f"{__class__.__name__} ERROR: 'files' list is empty")

        files.sort()

        if progress:
            files = tqdm.tqdm(files)

        actin = []
        for file in files:
            time.sleep(0.01)

            # Read spectrum:
            spec = self.ReadSpec(file, verb=verb, **spec_kw).spec

            # spec is the spectrograph class
            spectrum = spec.spectrum
            headers = spec.headers

            if not isinstance(spectrum, dict):
                continue

            if not indices:
                return pd.DataFrame([headers])

            if obj_in:
                headers['obj_in'] = obj_in
            
            # Calculate indices:
            calc_ind = self.CalcIndices(
                    spectrum, headers, table_df=table_df, indices=indices, verb=verb, **calcind_kw)
            

            data = dict()
            data.update(headers)
            data.update(calc_ind.indices)

            data['actin_ver'] = self.version

            actin.append(data)


        df = pd.DataFrame(actin)


        if save_data:
            if os.path.dirname(save_data):
                if not os.path.isdir(os.path.dirname(save_data)):
                    os.makedirs(os.path.dirname(save_data))
            df.to_csv(save_data, index=False)

        return df



    def plot_index_lines(self, file, index, table_df=None, wkey='wave', fkey='flux', ntol=2, show=True, spec_kw=dict(), **plt_kw):
        """Plot spectral lines and bandpasses used to calculate index 'index'. 
        """
        if isinstance(file, (list, np.ndarray, tuple)) and len(file) == 1:
            file = file[0]
        elif isinstance(file, (list, np.ndarray, tuple)) and len(file) > 1:
            raise RuntimeError(f"{sys._getframe().f_code.co_name} ERROR: Con only read one file at a time")


        if not os.path.isfile(file):
            raise FileNotFoundError(f"{sys._getframe().f_code.co_name} ERROR: File '{file}' not found")


        read_spec = self.ReadSpec(file, **spec_kw)
        
        if table_df is None:
            table = self.IndTable().table
        else:
            if not isinstance(table_df, pd.DataFrame):
                raise RuntimeError(f"{sys._getframe().f_code.co_name} ERROR: 'table_df' is not a pandas DataFrame")
            table = table_df

        if index not in np.unique(table.ind_id.values):
            raise RuntimeError(f"{sys._getframe().f_code.co_name} ERROR: Index {index} not available in 'ind_table'. Available indices are: {np.unique(table.ind_id.values)}")
        
        ind_tab = table[table.ind_id == index]
        ind_vars = ind_tab.ind_var.values


        _, axes = plt.subplots(1, len(ind_vars), figsize=(4*len(ind_vars), 4))

        for ind_var, ax in zip(ind_vars, axes):
            ctr = ind_tab[ind_tab.ind_var==ind_var].ln_ctr.values[0]
            win = ind_tab[ind_tab.ind_var==ind_var].ln_win.values[0]
            bt = ind_tab[ind_tab.ind_var==ind_var].bandtype.values[0]


            f = read_spec.spec.spectrum[fkey]
            w = read_spec.spec.spectrum[wkey]

            if len(f.shape) == 2:
                order = self.CalcIndices._spec_order(w, f, None, ctr, win, bt, type='dist', show_orders=False, verb=False)

                w = w[order]
                f = f[order]


            mask = (w > ctr - win * ntol)
            mask &= (w < ctr + win * ntol)
            
            f = f[mask]
            w = w[mask]

            ax.plot(w, f, **plt_kw)

            # add bandpasses:
            if bt == 'tri':
                bp = 1 - np.abs(w - ctr)/win
                bp = np.where(bp > 0, bp, bp*0.0)
            elif bt == 'sq':
                bp_mask = (w >= ctr - win/2) & (w <= ctr + win/2)
                bp = np.where(bp_mask, 1, 0.0)

            ax.plot(w, 1.2*f.max()*bp, 'k--', lw=0.7, label='bandpass')

            ax.axvline(ctr, color='k', ls=':', lw=0.7)
            ax.set_xlim(w.min(), w.max())
            #ax.set_ylim(f.min() - 0.2*f.min(), f.max() + 0.05*f.max())
            ax.set_xlabel("Wavelength [angstrom]")

        plt.tight_layout()
        if show:
            plt.show()



    def plot_index_line(self, file, line_id, table_df=None, wkey='wave', fkey='flux', show=True, ntol=4, show_bandpass=True, spec_kw=dict(), **plt_kw):
        """Plot spectral line 'line_id'.
        """
        if isinstance(file, (list, np.ndarray, tuple)) and len(file) == 1:
            file = file[0]
        elif isinstance(file, (list, np.ndarray, tuple)) and len(file) > 1:
            raise RuntimeError(f"{sys._getframe().f_code.co_name} ERROR: Con only read one file at a time")

        if not os.path.isfile(file):
            raise FileNotFoundError(f"{sys._getframe().f_code.co_name} ERROR: File '{file}' not found")
        
        read_spec = self.ReadSpec(file, **spec_kw)
        
        if table_df is None:
            table = self.IndTable().table
        else:
            if not isinstance(table_df, pd.DataFrame):
                raise RuntimeError(f"{sys._getframe().f_code.co_name} ERROR: 'table_df' is not a pandas DataFrame")
            table = table_df

        if line_id not in np.unique(table.ln_id.values):
            raise RuntimeError(f"{sys._getframe().f_code.co_name} ERROR: Line {line_id} not available in 'ind_table'. Available lines are: {np.unique(table.ln_id.values)}")

        ind_tab = table[table.ln_id == line_id]

        ctr = ind_tab.ln_ctr.values[0]
        win = ind_tab.ln_win.values[0]
        bt = ind_tab.bandtype.values[0]

        f = read_spec.spec.spectrum[fkey]
        w = read_spec.spec.spectrum[wkey]

        if len(f.shape) == 2:
            w, f, _ = self.CalcIndices._spec_order(w, f, None, ctr, win, bt, type='dist', show_orders=False, verb=False)

        mask = (w > ctr - win * ntol) 
        mask &= (w < ctr + win * ntol)
        
        f = f[mask]
        w = w[mask]

        plt.plot(w, f, **plt_kw)

        if show_bandpass:
            # add bandpasses:
            if bt == 'tri':
                bp = 1 - np.abs(w - ctr)/win
                bp = np.where(bp > 0, bp, bp*0.0)
            elif bt == 'sq':
                bp_mask = (w >= ctr - win/2) & (w <= ctr + win/2)
                bp = np.where(bp_mask, 1, 0.0)

            plt.plot(w, 1.2*f.max()*bp, 'k--', lw=0.7, label='bandpass')

            plt.axvline(ctr, color='k', ls=':', lw=0.7)
        plt.xlim(w.min(), w.max())
        plt.xlabel("Wavelength [angstrom]")

        plt.tight_layout()
        if show:
            plt.show()

        return ctr, win

    
