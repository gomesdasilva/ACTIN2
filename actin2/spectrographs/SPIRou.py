import sys, os
import numpy as np
import glob
import pandas as pd

from astropy.io import fits

from ._spec_tools import printif, read_fits, read_headers, wave_star_rest_frame, wave_corr_berv



# spec_hdrs_0 = dict(
#     obj     = 'OBJECT',
#     instr   = 'INSTRUME',
#     exptime = 'EXPTIME', # [sec]
#     gain = hdr['GAIN'], # [e-/ADU]
#     rdnoise = hdr['RDNOISE'], # [e-]
#     pi_name = hdr['PI_NAME'],
#     run_id = hdr['RUNID'],
#     date_obs = hdr['DATE-OBS'] + 'T' + hdr['UTIME'],
#     ra = hdr['RA_DEG'],
#     dec = hdr['DEC_DEG'],
#     airmass = hdr['AIRMASS'], # airmass at start
#     drs_ver = hdr['VERSION'],
# )

# spec_hdrs_1 = dict(
#     bjd = 'BJD',
#     berv = 'BERV', # [km/s]
#     snr25 = 'SNR25'
# )



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

    def __init__(self, file, hdu, obj_in=None, verb=False):
        printif("Running SPIRou", verb)

        spec, headers = self.read_spec_spirou(file, hdu)

        # correct barycentric motion (Etienne in email, s files)
        spec['wave_raw'] = wave_corr_berv(spec['wave_raw'], headers['berv'])

        #print(headers['rv'])

        # shift wave to target rest frame: (Etienne in email, s files)
        try:
            spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['rv'])
        except KeyError:
            try:
                spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['targ_rv'])
            except KeyError:
                printif("WARNING: wavelength not shifted to target rest frame", verb)

        #headers['noise'] = headers['gain'] * headers['rdnoise']

        spec['flux_err'] = np.sqrt(abs(spec['flux_raw']))
        # print(spec['flux_err'][40][:]) #! gives nan??

        # print(spec['flux_err'])

        # sys.exit()

        try:
            spec['flux'] = spec['flux_tell_corr']
        except KeyError:
            try:
                spec['flux'] = spec['flux_raw']/spec['blaze']
                #print("Flux was deblazed")
            except KeyError:
                spec['flux'] = spec['flux_raw']
                #print("Flux was NOT debalzed")



        self.spectrum = spec
        self.headers = headers

        #sys.exit()


    def read_spec_spirou(self, file, hdu):
        #! should check file type "s"?
        if hdu:
            pass
        else:
            hdu = fits.open(file)

        #print(file)
        #sys.exit()

        #print(hdu.info())
        #sys.exit()
        #print(hdu[0].header)
        #sys.exit()

        spec = dict()

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
        #print(hdu[1].header[:10])#;sys.exit()

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

            #snr25 = hdu[1].header['SNR25'], # snr at order 25

            berv = hdu[1].header["BERV"] * 1000,
            bjd = hdu[1].header['BJD'],
        )

        try:
            headers['targ_rv'] = hdr['OBJRV'] * 1e3 # [m/s]
        except KeyError:
            pass

        # get the CCF file:
        filename = hdr['FILENAME'][:-1] + "v.fits"

        file_v = os.path.join(os.path.dirname(file), filename)

        #* Make try statement:
        hdu = fits.open(file_v)
        #print(hdu.info())
        #print(hdu['CCF'].header) #! CCF data
        #print(hdu['CCF'].data.names)

        ccf_rv = np.zeros(len(hdu['CCF'].data))
        ccf_profile = np.zeros(len(hdu['CCF'].data))
        for i, data in enumerate(hdu['CCF'].data):
            ccf_rv[i] = data['Velocity']
            ccf_profile[i] = data['combined']

        spec['ccf_rv'] = ccf_rv
        spec['ccf_profile'] = ccf_profile

        #headers['rv'] = hdu[0].header['CCFRVC'] * 1e3

        try:
            headers['rv'] = hdu[0].header['CCFRVC'] * 1e3
            headers['rv_phn'] = hdu[0].header['DVRMS']
            headers['fwhm'] = hdu[0].header['CCFFWHM'] * 1e3
            headers['cont'] = hdu[0].header['CCFCONT']
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





        
