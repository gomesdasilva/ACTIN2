import sys, os
import numpy as np
import glob
import pandas as pd

from astropy.io import fits

from ._spec_tools import read_hdr_data
from ._spec_tools import printif
from ._spec_tools import wave_star_rest_frame


SP_HDRS = dict(
    instr = 'INSTRUME',
    obj = 'OBJECT',
    bjd = 'HIERARCH ESO QC BJD',
    date_obs = 'DATE-OBS',
    exptime = 'EXPTIME',
    snr60 = 'HIERARCH ESO QC ORDER85 SNR',
    airmass_start = 'HIERARCH ESO TEL3 AIRM START',
    airmass_end = 'HIERARCH ESO TEL3 AIRM END',
    seeing_start = 'HIERARCH ESO TEL3 AMBI FWHM START', # [arcsec]
    seeing_end = 'HIERARCH ESO TEL3 AMBI FWHM END', # [arcsec]
    ra = 'RA',      # [deg]
    dec = 'DEC',     # [deg]
    prog_id = 'HIERARCH ESO OBS PROG ID',
    pi_coi = 'PI-COI',
    ftype = 'HIERARCH ESO PRO CATG',
    drs_id = 'HIERARCH ESO PRO REC1 DRS ID',
    targ_rv = 'HIERARCH ESO OCS OBJ RV',
    rv_step = 'HIERARCH ESO RV STEP', #=     0.5 / Coordinate increment per pixel
    rv_ref_pix = 'HIERARCH ESO RV START', #= -21.11 / Coordinate at reference pixel
    berv = 'HIERARCH ESO QC BERV', # [km/s]
    rv = 'HIERARCH ESO QC CCF RV', #= -1.19840358654132 / Radial velocity [km/s]        
    rv_err = 'HIERARCH ESO QC CCF RV ERROR', #= 0.00053474023623152 / Uncertainty on radial veloc
    fwhm = 'HIERARCH ESO QC CCF FWHM', #= 9.4897932752335 / CCF FWHM [km/s]                  
    fwhm_err = 'HIERARCH ESO QC CCF FWHM ERROR', #= 0.00106948047246304 / Uncertainty on CCF FWHM [
    cont = 'HIERARCH ESO QC CCF CONTRAST', #= 59.8806426078076 / CCF contrast %                
    cont_err = 'HIERARCH ESO QC CCF CONTRAST ERROR', #= 0.00674842708267666 / CCF contrast error % 
    continuum = 'HIERARCH ESO QC CCF CONTINUUM', #= 2560440.01971191 / CCF continuum level [e-]     
    mask = 'HIERARCH ESO QC CCF MASK', #= 'F9      ' / CCF mask used                           
    flux_asym = 'HIERARCH ESO QC CCF FLUX ASYMMETRY', #= 0.0126966418595723 / CCF asymmetry (km/s)  
    flux_asym_err = 'HIERARCH ESO QC CCF FLUX ASYMMETRY ERROR', #= 0.000921109029304877 / CCF asymmetry 
    bis = 'HIERARCH ESO QC CCF BIS SPAN', #= -0.01440933698031 / CCF bisector span (km/s)     
    bis_err = 'HIERARCH ESO QC CCF BIS SPAN ERROR', #= 0.00106948047246304 / CCF bisector span err
    blaze_file = "HIERARCH ESO PRO REC2 CAL12 NAME"
)




class ESPRESSO:

    def __init__(self, file, hdu, obj_in, rv_in, verb, **spec_kw):

        printif("Running ESPRESSO", verb)
        
        spec, hdr = self.read_fits_espresso(hdu=hdu)

        headers = read_hdr_data(hdr, SP_HDRS)

        for h in headers:
            if h in ['rv', 'fwhm', 'bis', 'berv']:
                headers[h] = headers[h]*1000

        # shift wave to target rest frame:
        spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['rv'])

        # TODO: add CCF profile


        # import matplotlib.pylab as plt
        # if 'S2D' in headers['ftype']:
        #     plt.plot(spec['wave'][12], spec['flux_raw'][12], label="wave")
        #     #plt.plot(spec['wave'][15], spec['flux_raw'][15], label="wave")
        # if 'S1D' in headers['ftype']:
        #     plt.plot(spec['wave'], spec['flux_raw'], label="flux_raw", marker='.')
        #     #plt.plot(spec['wave'], spec['flux_cal'], label="flux_cal")
        #     #plt.plot(spec['wave'], spec['flux_cal_skysub'], label='flux_cal_skysub')
        #     #plt.plot(spec['wave_raw'], spec['flux_raw'], label="wave_raw")
        # plt.legend(loc=1)

        # from actin2 import IndTable
        # df = IndTable().table
        # ctr = df[df.ln_id == 'CaIIH'].ln_ctr.values[0]
        # win = df[df.ln_id == 'CaIIH'].ln_win.values[0]
        # plt.axvline(df[df.ln_id == 'CaIIH'].ln_ctr.values[0], c='k')
        # plt.axvline(ctr - win, c='k', ls='--')
        # plt.axvline(ctr + win, c='k', ls='--')

        # plt.show()

        # #print(hdr)
        # sys.exit()


        spec['flux'] = spec['flux_raw']
        headers['spec_flg'] = 'OK'

        max_RON = 8 # [e-/pix] from the pipeline manual, pag. 16
        CONAD = 1.1 # [e-/ADU] from the pipeline manual, pag. 16
        headers['ron'] = max_RON * CONAD

        self.spectrum = spec
        self.headers = headers
        self._headers = headers

    def read_fits_espresso(self, file=None, hdu=None):

        if not hdu:
            hdu = fits.open(file)

        hdr = hdu[0].header

        if hdr['HIERARCH ESO PRO CATG'] == 'S1D_A':
            spec = dict(
                flux_raw = hdu[1].data['flux'], # debalzed
                flux_err = hdu[1].data['error'],
                wave_raw = hdu[1].data['wavelength_air'], # already in target rest frame
        )

        elif hdr['HIERARCH ESO PRO CATG'] == 'S2D_A':
            spec = dict(
                flux_raw = hdu[1].data, # blazed
                flux_err = hdu[2].data,
                wave_raw = hdu[5].data # wavelength air (bary)
            )
        
        else:
            raise ValueError(f"ERROR: File '{hdr['HIERARCH ESO PRO CATG']}' not implemented")

        # if hdr['HIERARCH ESO PRO CATG'] == 'S1D_FINAL_A':
        #     #! The errors are not good!
        #     # notes: flux_el is zero
        #     spec = dict(
        #         flux_raw = hdu[1].data['flux'][0], # blazed
        #         flux_err = hdu[1].data['err'][0],
        #         flux_cal = hdu[1].data['flux_cal'][0],
        #         flux_cal_skysub = hdu[1].data['FLUX_CAL_SKYSUB'][0],
        #         wave_raw = hdu[1].data['wave_air'][0] # wavelength air (bary)
        #     )

        hdu.close()

        return spec, hdr

