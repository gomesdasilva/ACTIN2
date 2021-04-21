import sys, os
import numpy as np
import glob
import pandas as pd

from astropy.io import fits

from ._spec_tools import read_hdr_data
from ._spec_tools import printif
from ._spec_tools import wave_star_rest_frame
from ._spec_tools import wave_corr_berv




SP_HDRS = dict(

)


class SPIRou:
    """
    Notes:
    ------
    - HDU has 49 orders, 4088 pixels per order.
    - Wavelength resolution is not constant if using wave data.
    - hdu[0].header['CCFRVC'] only in "v" files (CCF).

    - "s" (science) files are 1d deblazed
    - hdu[1] is constant wave resolution of 0.05 ang
    """

    # TODO: Include 2D spectra

    def __init__(self, file, hdu, obj_in, rv_in, verb, **spec_kw):
        printif("Running SPIRou", verb)

        spec, headers = self.read_spec_spirou(file, hdu)

        # correct barycentric motion (Etienne in email)
        spec['wave_raw'] = wave_corr_berv(spec['wave_raw'], headers['berv'])

        # shift wave to target rest frame: (Etienne in email)
        try:
            spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['rv'])
        except KeyError:
            try:
                spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['targ_rv'])
            except KeyError:
                print("WARNING: wavelength not shifted to target rest frame")

        headers['ron'] = headers['gain'] * headers['rdnoise']

        try:
            spec['flux'] = spec['flux_tell_corr']
        except KeyError:
            try:
                print("Flux was deblazed")
                spec['flux'] = spec['flux_wave']/spec['blaze']
            except KeyError:
                print("Flux was NOT debalzed")
                spec['flux'] = spec['flux_raw']



        self.spectrum = spec
        self.headers = headers
        self._headers = headers

        #sys.exit()


    def read_spec_spirou(self, file, hdu):
        #! should check file type "s"?
        if hdu:
            pass
        else:
            hdu = fits.open(file)

        print(file)
        #sys.exit()

        print(hdu.info())
        #sys.exit()
        print(hdu[0].header)
        #sys.exit()

        if hdu[1].header['EXTNAME'] == 'UniformWavelength':
            spec = dict(
                wave_raw = hdu[1].data.field(0) * 10, # to ang
                flux_raw = hdu[1].data.field(1),
                flux_tell_corr = hdu[1].data.field(9)
            )

        if hdu[1].header['EXTNAME'] == 'FluxAB':
            spec = dict(
                flux_raw = hdu[1].data[:],
                wave_raw = hdu[2].data[:] * 10,
                blaze = hdu[3].data[:],
                flux_atm = hdu[4].data[:]
            )

        hdr = hdu[0].header

        headers = dict(
            obj = hdr['OBJECT'],
            instr = hdr['INSTRUME'],
            exptime = hdr['EXPTIME'], # [sec]
            gain = hdr['GAIN'], # [e-/ADU]
            rdnoise = hdr['RDNOISE'], # [e-]
            pi_name = hdr['PI_NAME'],
            run_id = hdr['RUNID'],
            date_obs = hdr['DATE-OBS'] + 'T' + hdr['UTIME'],
            ra = hdr['RA_DEG'],
            dec = hdr['DEC_DEG'],
            airmass = hdr['AIRMASS'], # airmass at start
            drs_ver = hdr['VERSION'],

            snr25 = hdu[1].header['SNR25'], # snr at order 25

            berv = hdu[1].header["BERV"] * 1000,
            bjd = hdu[1].header['BJD'],
        )

        try:
            headers['targ_rv'] = hdr['OBJRV'] * 1000 # [m/s]
        except KeyError:
            pass

        filename = hdr['FILENAME'][:-1] + "v.fits"

        file_v = os.path.join(os.path.dirname(file), filename)

        hdu = fits.open(file_v)

        try:
            headers['rv'] = hdu[1].header['CCFRVC'] * 1000
        except KeyError:
            print("ERROR: No 'CCFRVC' header on hdu[1]")



        #file_e:
        #EXTNAME == 'FluxAB'
        # hdu = fits.open('1234567e.fits')
        # # looking at science data channel (AB)
        # snr_order34 = hdu[1].header['SNR34']
        # flux_order34 = hdu[1].data[34]
        # wavelength_order34 = hdu[5].data[34]

        #file_v = file
        #EXTNAME == "CCF"

        #file_s:
        #EXTNAME == 'UniformWavelength' 

        #rv = hdr['CCFRVC'] * 1000,


        hdu.close()

        return spec, headers





        
