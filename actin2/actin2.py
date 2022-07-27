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

#? important: when using 1d files, the flux is already subtracted by the blaze function and therefore the sqrt(F) photon errors lose the physical meaning.

# TODO: add option/class/function to plot lines from IntTable (with the bandpass)

# TODO POST-PROCESSING: option to select HARPS-pre/HARPS-pos (option)
# TODO POST-PROCESSING: add correction for CCF FWHM and contrast drift (option)?

class ACTIN:
    """The ACTIN class. Reads fits files and calculates activity indices.

    Attributes:
        ReadSpec (actin2.ReadSpec) : Object that reads spectrum and   
            headers.
        IndTable (actin2.IndTable) : Object containing the 
            indices parameters table.
        CalcIndices (actin2.CalcIndices) : Object that calculates the 
            indices.
    """

    def __init__(self):
        self.version = '2.0.0 beta 5'
        self.ReadSpec = ReadSpec
        self.IndTable = IndTable
        self.CalcIndices = CalcIndices

    #@timeit
    def run(self, files, indices, obj_in=None, save_data=None, progress=False, verb=False, table_df=None, **spec_kw):
        """Run ACTIN for a list of files.

        Args:
            files (list) : List of fits files paths to be read.
            indices (list) : List of indices identification (as in 
                the indices table) to be calculated. If empty returns the fits headers only.
            rv_in (float, None) : If float, use this RV value to 
                shift spectrum to stellar rest frame, if None use RV value from fits file.
            obj_in (str, None) : If not None, the target name from 
                the fits file is overriden by this input.
            save_data (str, None): path with file_name to save output in csv format.
            progress (bool) : If True show progress bar.
            table_df (None, pd.DataFrame) : Table with the indices parameters.
                If None use the built-in table.
            verb (bool) : Activate verbose option.
        
        Returns:
            pandas.DataFrame : Output pandas table.
        """

        printif("\nRUNNING ACTIN 2", verb)
        printif("---------------", verb)

        # Enforce 'files' type as list:
        if not isinstance(files, (list, np.ndarray, tuple)):
            files = [files]

        assert len(files) != 0, "'files' list is empty"

        files.sort()

        if progress:
            verb = False
            files = tqdm.tqdm(files)

        actin = [] 
        for file in files:
            time.sleep(0.01)
            read_spec = self.ReadSpec(file, verb=verb, **spec_kw)

            if not read_spec.spec.spectrum:
                continue

            if not indices:
                return pd.DataFrame([read_spec.spec.headers])

            # read_spec.spec is the spectrograph class
            # TODO: change inputs below to wave, flux, etc?
            try:
                calc_ind = self.CalcIndices(read_spec.spec, table_df=table_df, indices=indices, interp=True, verb=verb)
            except:
                print("ERROR: calc_ind")
                continue

            data = calc_ind.data

            # TODO: make function to read version from VERSION
            data['actin_ver'] = self.version

            actin.append(data)


        df = pd.DataFrame(actin)

        # TODO: option to save rdb files?

        if save_data:
            if not os.path.isdir(os.path.dirname(save_data)):
                os.makedirs(os.path.dirname(save_data))
            df.to_csv(save_data, index=False)
            printif(f"Data saved to: {save_data}", verb)

        return df



    # TODO: Test 2D spectrum
    def plot_index_lines(self, file, index, table_df=None, wkey='wave', fkey='flux', show=True):
        """Plot sepctral lines and bandpasses used to calculate index 'index'. 
        """
        read_spec = self.ReadSpec(file)
        
        if table_df is None:
            table = self.IndTable().table
        else:
            table = table_df
        
        ind_tab = table[table.ind_id == index]
        ind_vars = ind_tab.ind_var.values

        fig, axes = plt.subplots(nrows=1, ncols=len(ind_vars), figsize=(4*len(ind_vars), 4))

        for ind_var, ax in zip(ind_vars, axes):
            ctr = ind_tab[ind_tab.ind_var==ind_var].ln_ctr.values[0]
            win = ind_tab[ind_tab.ind_var==ind_var].ln_win.values[0]
            bt = ind_tab[ind_tab.ind_var==ind_var].bandtype.values[0]

            if bt == 'sq': win /= 2

            f = read_spec.spec.spectrum[fkey]
            w = read_spec.spec.spectrum[wkey]

            mask = (w > ctr - win - 2 * win) 
            mask &= (w < ctr + win + 2 * win)
            
            f = f[mask]
            w = w[mask]

            ax.plot(w, f/f.mean())

            ax.axvline(ctr, color='k', ls=':', lw=0.7)
            ax.axvline(ctr - win, color='k', ls='--', lw=0.7)
            ax.axvline(ctr + win, color='k', ls='--', lw=0.7)
            ax.set_xlim(w.min(), w.max())
            ax.set_xlabel("Wavelength [angstrom]")

        plt.tight_layout()
        if show:
            plt.show()



    def plot_index_line(self, file, line_id, table_df=None, wkey='wave', fkey='flux', show=True, ntol=2, show_bandpass=True, **kwargs):
        """Plot sepctral line 'line_id'.
        """
        read_spec = self.ReadSpec(file)
        
        if table_df is None:
            table = self.IndTable().table
        else:
            table = table_df

        ind_tab = table[table.ln_id == line_id]

        ctr = ind_tab.ln_ctr.values[0]
        win = ind_tab.ln_win.values[0]
        bt = ind_tab.bandtype.values[0]

        if bt == 'sq': win /= 2

        f = read_spec.spec.spectrum[fkey]
        w = read_spec.spec.spectrum[wkey]

        mask = (w > ctr - win - win * ntol) 
        mask &= (w < ctr + win + win * ntol)
        
        f = f[mask]
        w = w[mask]

        plt.plot(w, f/f.mean(), **kwargs)

        if show_bandpass:
            plt.axvline(ctr, color='k', ls=':', lw=0.7)
            plt.axvline(ctr - win, color='k', ls='--', lw=0.7)
            plt.axvline(ctr + win, color='k', ls='--', lw=0.7)
        #plt.xlim(ctr - win - win*2, ctr + win + win*2)
        plt.xlim(w.min(), w.max())
        plt.xlabel("Wavelength [angstrom]")

        plt.tight_layout()
        if show:
            plt.show()

        return ctr, win

    
