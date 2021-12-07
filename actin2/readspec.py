import os
import importlib
import matplotlib.pylab as plt
from astropy.io import fits

from . import spectrographs
from .spectrographs._spec_tools import printif




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
    # TODO 6: Add SPIRou
    # OPTIONAL:
    # TODO: Include barycorr.py!
    # TODO: add option to input RV and BERV!
    # TODO: Add option to retrieve spectra without being at rest frame!

    def __init__(self, file, obj_in=None, verb=False, spec_class_in=None, **spec_kw):
        printif("Running ReadSpec", verb)

        # check if file exists:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"*** ERROR: File not found")


        printif(f"loaded file: {file}", verb)

        # Get instrument:
        hdu = fits.open(file)
        try:
            instr = hdu[0].header['INSTRUME']
        except KeyError:
            try:
                instr = hdu[0].header['HIERARCH ESO OBS INSTRUMENT']
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
                spectrograph = importlib.import_module("." + instr, "actin2.spectrographs").__getattribute__(instr)
            except ModuleNotFoundError:
                raise ModuleNotFoundError(f"*** ERROR: Instrument '{instr}' not detected in 'spectrographs'. Available instruments are: {spectrographs.__all__}")


        # Create spectrograph object, this object will have the attributes given to him by the spectrograph class (can be different for different spectrographs: e.g. having bisector dictionary):
        self.spec = spectrograph(file, hdu, obj_in, verb, **spec_kw)

        self.spec.headers['file'] = os.path.basename(file)

        hdu.close()

        #* Required attributes:
        #self.spectrum = spectrograph.spectrum
        #self.headers = spectrograph.headers   # headers for output



    def plot_spec(self, key_wave='wave', key_flux='flux', order=None, ax=None, show=False, **plt_kw):
        """Plot spectrum.

        Args:
            key_wave (str, option) : Keyword for the wavelength array stored 
                in ``spectrum`` dictionary.
            key_flux (str, option) : Keyword for the flux array stored 
                in ``spectrum`` dictionary.
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


