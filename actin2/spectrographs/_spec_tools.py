"""
ACTIN 2: General functions used to process data from spectrographs.
"""
import numpy as np
import pandas as pd
from astropy.io import fits

# THIS IS A GENERAL FUNCTION, SHOULD GO TO ANOTHER FILE?
def printif(text, verb=True):
    """Print text if verb is True"""
    if verb:
        print(text)


def wave_star_rest_frame(wave, rv):
    """Change wavelength to stellar rest frame
    'wave' must be in the solar system baricentric frame."""
    c = 299792458.0 # [m/s]
    dwave = rv * wave / c
    wave_corr = wave - dwave
    return wave_corr


def wave_corr_berv(wave, berv):
    """Change wavelength to solar system baricentric frame"""
    c = 299792458.0 # [m/s]
    dwave = - berv * wave / c
    wave_corr = wave - dwave
    return wave_corr


def filter_headers(OUT_HDR_KEYS, all_hdr_dict):
    headers = dict()
    for key in OUT_HDR_KEYS:
            if key in all_hdr_dict:
                headers[key] = all_hdr_dict[key]
    return headers


def read_fits(fits_file=None, hdu=None, instr=None, calc_x=True):
    if fits_file:
        hdu = fits.open(fits_file)
    y = hdu[0].data
    hdr = hdu[0].header

    hdu.close()

    if calc_x:
        if instr in ['HARPS', 'CORALIE']:
            obs = 'ESO'
        if instr == 'HARPN':
            obs = 'TNG'
        try:
            x = calc_fits_x_2d(hdr, obs)
        except KeyError:
            try:
                x = calc_fits_x_1d(hdr)
            except:
                print("*** Error reading fits x coordinate")
                x = np.nan
        return x, y, hdr
    else:
        return y, hdr


def calc_fits_x_1d(hdr, key_a='CRVAL1', key_b='CDELT1', key_c='NAXIS1'):
    return hdr[key_a] + hdr[key_b] * np.arange(hdr[key_c])


def calc_fits_x_2d(hdr, obs):
    deg = hdr[f'HIERARCH {obs} DRS CAL TH DEG LL'] # polynomial degree

    ll_coeff = np.zeros((hdr['NAXIS2'], deg + 1))

    # Read coefficients
    for i in range(hdr['NAXIS2']):
        for j in range(deg + 1):
            ll_coeff[i, j] = hdr['HIERARCH {} DRS CAL TH COEFF LL{}'.format(obs, (j + (deg + 1)*i))]

    # Evaluate polynomials
    x = np.arange(hdr['NAXIS1'])  # Pixel array
    wave_raw = np.zeros([hdr['NAXIS2'], hdr['NAXIS1']])  # Wavelength 2D array
    for i in range(len(wave_raw)):
        wave_raw[i] = np.poly1d(ll_coeff[i][::-1])(x)
    return wave_raw




def read_headers(hdr, headers, data=None, verb=False):
    """Read fits header data using 'headers' dictionary.
    Result is included in 'data' dictionary if not 'None'. If 'None' a new dictionary is returned"""
    if not data:
        data = {}

    for key, hdr_id in zip(headers.keys(), headers.values()):
        try:
            data[key] = hdr[hdr_id]
        except KeyError:
            printif(f"Header {key} not in fits file", verb)
            data[key] = None
    return data


# 2.0_beta_12
import astropy.units as u
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler


def make_specutils(wave, flux, flux_e):
    from specutils import Spectrum1D
    from astropy.nddata import StdDevUncertainty
    import astropy.units as u
    if flux_e is not None:
        return Spectrum1D(spectral_axis=wave*u.AA, flux=flux*u.dimensionless_unscaled, uncertainty=StdDevUncertainty(flux_e*u.dimensionless_unscaled))
    else:
        return Spectrum1D(spectral_axis=wave*u.AA, flux=flux*u.dimensionless_unscaled)


def from_specutils(spec):
    try:
        return spec.wavelength.value, spec.flux.value, spec.uncertainty.array
    except AttributeError:
        return spec.wavelength.value, spec.flux.value
    


def interpolate_region(wave: np.array, flux: np.array, flux_e: np.array, wmin: float, wmax: float, grid='arange', n_points=10, wstep=0.01, method='cubic'):
    """Interpolate spectral region.

    Parameters:
    -----------
        wave : np.array
            Spectral wavelength (A).
        flux : np.array
            Spectral flux.
        flux_e : np.array
            Spectral flux uncertainty.
        wmin : float
            Lower limit of wavelength region (A).
        wmax : float
            Upper limit of wavelength region (A).
        grid : str
            Type of wavelength grid: if 'arange' (default) use numpy.arange function with a step of 'wstep'; if 'linspace' use numpy.linspace function with 'n_points' number of points.
        n_points : int
            Number of points in the interpolated region if using grid type 'linspace'.
        wstep : float
            Wavelength step used in 'arange' grid (A).
        method : str
            Method used for interpolation: if 'cubic' (default) use cubic interpolation; if 'linear' use linear interpolation.
    
    Returns
    -------
        wave_re : np.array
            Resampled wavelwngth in region (A)
        flux_re : np.array
            Resampled flux in region (arbitray units)
        flux_e_re : np.array
            Resampled flux uncertainty in region (arbitrary units)
    """
    if grid=='linspace':
        x_interp = np.linspace(wmin, wmax, n_points)
    if grid=='arange':
        x_interp = np.arange(wmin, wmax, wstep)
    if flux_e is not None:
        mask = np.isnan(flux) | np.isnan(flux_e)
        wave = wave[~mask]
        flux = flux[~mask]
        flux_e = flux_e[~mask]
        spec = make_specutils(wave, flux, flux_e)
    else:
        mask = np.isnan(flux)
        wave = wave[~mask]
        flux = flux[~mask]
        spec = make_specutils(wave, flux, None)
    if method == 'cubic':
        resampler = SplineInterpolatedResampler(bin_edges='nan_fill')
    if method == 'linear':
        resampler = LinearInterpolatedResampler(extrapolation_treatment='nan_fill')
    spec_re = resampler(spec, x_interp * u.AA)
    return from_specutils(spec_re)