import sys, os
import importlib
from importlib.machinery import SourceFileLoader
import matplotlib.pylab as plt
from astropy.io import fits
import numpy as np

from . import spectrographs
from .spectrographs._spec_tools import printif



class ReadSpec:
    """
    Extracts data from spectrograph fits files.

    Args:
        file (str): Fits file path containing the spectral data.
        spec_class_in (str, None): Spectrograph class identification. Useful if new spectrographs are added.
        verb (bool): Turn verbose on/off.
        **spec_kw (dict): Additional keyword arguments to be passed to the spectrograph class.

    *Class attributes:*

    Attributes:
        spectrum (dict) : Dictionary containing wavelength and flux
        header (dict) : Dictionary containing selected headers from 
            the fits file
        spec (object): The spectrograph class.

    """

    def __init__(self, file, verb=False, spec_class_in=None, spec_file_in=None, **spec_kw):

        # check if file exists:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"{__class__.__name__} ERROR: File not found")


        #printif(f"loaded file: {file}", verb)

        hdu = fits.open(file)

        # using a spectrograph file outside ACTIN path:
        if spec_file_in:
            instr = os.path.basename(spec_file_in).split('.')[0]
            spectrograph = SourceFileLoader(instr, spec_file_in).load_module().__getattribute__(instr)

        else:
            # the spectrograph class is given as an input (to add new spectrographs easily):
            if spec_class_in:
                instr = spec_class_in
            else:
                try:
                    instr = hdu[0].header['INSTRUME']
                except KeyError:
                    try:
                        instr = hdu[0].header['HIERARCH ESO OBS INSTRUMENT']
                    except KeyError:
                        for instrument in spectrographs.__all__:
                            if instrument in os.path.basename(file):
                                instr = instrument

            #printif(f"loaded instr: {instr}", verb)


            # import class from module with names 'instr': the instrument class
            try:
                spectrograph = importlib.import_module("." + instr, "actin2.spectrographs").__getattribute__(instr)
            except ModuleNotFoundError:
                raise ModuleNotFoundError(f"{__class__.__name__} ERROR: Instrument '{instr}' not detected in 'spectrographs' directory. Available instruments are: {spectrographs.__all__}")

            #print(spectrograph);sys.exit()


        # Create spectrograph object, this object will have the attributes given to him by the spectrograph class (can be different for different spectrographs: e.g. having bisector dictionary):
        self.spec = spectrograph(hdu, file=file, verb=verb, **spec_kw)

        # Spectrum dictionary:
        self.spectrum = self.spec.spectrum
        # Headers dictionary:
        self.headers = self.spec.headers


        self.spec.headers['file'] = os.path.basename(file)

        hdu.close()



    def plot_spec(self, key_wave='wave', key_flux='flux', order=None, ax=None, show=False, **plt_kw):
        """Plot spectrum.

        Args:
            key_wave (str, option): Keyword for the wavelength array stored 
                in ``spectrum`` dictionary.
            key_flux (str, option): Keyword for the flux array stored 
                in ``spectrum`` dictionary.
            order (int, None): Spectral order to plot if using 2D spectrum.
            ax (matplotlib.axes, None): Axes for the plot. If ``None`` creates ``plt`` axes.
            show (bool): If ``True`` show the plot.
            **plt_kw: Additional keyword arguments to be passed to ``plt.plot()``.
        """
        spec = self.spec.spectrum

        if not ax:
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
            raise ValueError(f"{sys._getframe().f_code.co_name} ERROR: This is a 2D spectrum, 'order' value should be between 0 and {spec[key_flux].shape[0]-1} in the plot function")

        ax.plot(wave, flux, **plt_kw)

        if show:
            plt.show()

    
    def plot_ccf(self, order=-1, show=False, ax=plt, **plt_kw):

        try: 
            profile = self.spec.ccf_profile['profile'][order]
            rv = self.spec.ccf_profile['rv'] # [km/s]
        except AttributeError:
            raise AttributeError(f"{sys._getframe().f_code.co_name} ERROR: This spectrum has no CCF profile")

        ax.plot(rv, profile, **plt_kw)

        if show:
            plt.show()

    
    def plot_bisector(self, show=False, ax=plt, **plt_kw):

        try: 
            ccf_bisector = self.spec.ccf_bisector['bisector']
            ccf_rv = self.spec.ccf_bisector['rv'] # [km/s]
        except AttributeError:
            raise AttributeError(f"{sys._getframe().f_code.co_name} ERROR: This spectrum has no CCF bisector")

        ax.plot(ccf_rv, ccf_bisector, **plt_kw)

        if show:
            plt.show()


