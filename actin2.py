import sys, os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd

#* TEST CLASS:
#*--------------------
#* MAKE A CLASS WITH THESE FUNCTIONS WITH UPGRADABILITY IN MIND:
# - Ability to add, delete spectral lines and load, save different tables [done]
# - Ability to add new spectrographs easily
#    - Ex. each spectrograph can be included in a file with spectrograph name
# - Ability to change the function that calculates the flux easily




class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


import spectrographs

from spectrographs._spec_tools import printif

import importlib

import decorators

#from spectrographs import *


# read fits file and return spectrum + selected headers
class ReadSpec:
    """
    Extracts data from spectrograph fits files.

    Attributes:
        spectrum (dict) : Dictionary containing wavelength and flux
        header (dict) : Dictionary containing selected headers from 
            the fits file

    New features:

    - RON added to flux noise
    - (BERV correction with barycorr)
    - Extracts CCF profile (and Bisector)
    """
    # REQUIRED:
    # TODO: Add option to save spectrum, CCF and BIS - connect to ACTIN class
    # TODO 2: Add function to plot CCF Bisector
    # TODO 4: Add ESPRESSO 1D and 2D
    # TODO 5: Read ASCII table (wave, flux, err) -> make new file/class in spectrographs
    # TODO 6: Add SPIRou
    # OPTIONAL:
    # TODO: Include barycorr.py!

    def __init__(self, file, obj_in=None, verb=True, spec_class_in=None, **spec_kw):
        printif("Running ReadSpec", verb)

        # check if file exists:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"{bcolors.FAIL}READSPEC ERROR:{bcolors.ENDC} File not found")


        printif(f"loaded file: {file}", verb)

        # Get instrument:
        from astropy.io import fits
        hdu = fits.open(file)
        try:
            instr = hdu[0].header['INSTRUME']
        except KeyError:
            for instrument in spectrographs.__all__:
                if instrument in os.path.basename(file):
                    instr = instrument

        printif(f"loaded instr: {instr}", verb)

        # the spectrograph class is given as an input (to add new spectrographs easily):
        if spec_class_in:
            spectrograph = spec_class_in
        else:
            # import class from module with names 'instr': the instrument class
            try:
                spectrograph = importlib.import_module(f"spectrographs.{instr}").__getattribute__(instr)
            except ModuleNotFoundError:
                raise ModuleNotFoundError(f"{bcolors.FAIL}READSPEC ERROR:{bcolors.ENDC} Instrument '{instr}' not detected in 'spectrographs'. Available instruments are: {spectrographs.__all__}")


        # Create spectrograph object, this object will have the attributes given to him by the spectrograph class (can be different for different spectrographs: e.g. having bisector dictionary):
        self.spec = spectrograph(file, hdu, obj_in, verb, **spec_kw)

        self.spec.headers['file'] = os.path.basename(file)

        hdu.close()

        #* Required attributes:
        #self.spectrum = spectrograph.spectrum
        #self.headers = spectrograph.headers   # headers for output



    def plot(self, key_wave='wave', key_flux='flux', order=None, ax=None, show=False, **plt_kw):
        """Plot spectrum.

        Args:
            key_wave (str, option) : Keyword for the wavelength array stored 
                in ``spectrum`` dictionary.
            key_flux (str, option) : Keyword for the flux array stored 
                in ``spectrum`` dictionary.
        """
        spec = self.spec.spectrum

        if not ax:
            import matplotlib.pylab as plt
            ax = plt

        if order and len(spec[key_flux].shape) == 2:
            try:
                wave = spec[key_wave][order]
                flux = spec[key_flux][order]
            except IndexError as err:
                raise IndexError(err)

        elif len(spec[key_flux].shape) == 1:
            wave = spec[key_wave]
            flux = spec[key_flux]
        else:
            raise ValueError(f"This is a 2D spectrum, 'order' value should be between 0 and {spec[key_flux].shape[0]-1} in the plot function")

        ax.plot(wave, flux, **plt_kw)

        if show:
            plt.show()

    
    def plot_ccf(self, order=-1, show=False, ax=plt, **plt_kw):

        try: 
            profile = self.spec.ccf_profile['profile'][order]
            rv = self.spec.ccf_profile['rv'] # [km/s]
        except AttributeError:
            raise AttributeError("This spectrum has no CCF profile")

        ax.plot(rv, profile, **plt_kw)

        if show:
            plt.show()

    
    def plot_bisector(self, show=False, ax=plt, **plt_kw):

        try: 
            ccf_bisector = self.spec.ccf_bisector['bisector']
            ccf_rv = self.spec.ccf_bisector['rv'] # [km/s]
        except AttributeError:
            raise AttributeError("This spectrum has no CCF bisector")

        ax.plot(ccf_rv, ccf_bisector, **plt_kw)

        if show:
            plt.show()





#####


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

        self.table = pd.read_csv(table_csv)
        
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

####


####

import tqdm
import time

from CalcIndex import CalcIndices

from decorators import timeit

# TODO POST-PROCESSING: option to select HARPS-pre/HARPS-pos (option)
# TODO POST-PROCESSING: add correction for CCF FWHM and contrast drift (option)?

class ACTIN:
    """The ACTIN class. Reads fits files and calculates activity indices.

    Attributes:
        ReadSpec (actin2.ReadSpec) : Object that reads spectrum and   
            headers.
        IndTable (actin2.IndTable) : Object containing the 
            indices table.
        ProcessSpec (actin2.ProcessSpec) : Object to process the 
            spectrum.
    """

    def __init__(self):
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
            progress (bool) : If True show progress bar.
            table_df (None, pd.DataFrame) : Table with the indices parameters. If None use the built-in table.
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

            if not indices:
                return pd.DataFrame([read_spec.spec.headers])

            # read_spec.spec is the spectrograph class
            calc_ind = self.CalcIndices(read_spec.spec, table_df=table_df, indices=indices, interp=True, verb=verb)

            data = calc_ind.data

            data['actin_ver'] = '2.0.0'

            actin.append(data)


        df = pd.DataFrame(actin)

        if save_data:
            df.to_csv(save_data, index=False)
            printif(f"Data saved to: {save_data}", verb)

        return df

