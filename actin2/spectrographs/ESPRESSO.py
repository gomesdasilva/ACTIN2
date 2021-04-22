import sys, os
import numpy as np
import glob
import pandas as pd
from astropy.io import fits

from ._spec_tools import printif, read_fits, read_headers, wave_star_rest_frame, wave_corr_berv



spec_hdrs = dict(
    obj      = 'OBJECT',
    instr    = 'INSTRUME',
    date_obs = 'DATE-OBS',
    bjd      = 'HIERARCH ESO QC BJD',
    rv       = 'HIERARCH ESO QC CCF RV', # Radial velocity [km/s]
    rv_err   = 'HIERARCH ESO QC CCF RV ERROR', # Uncertainty on radial velocity [m/s]
    berv     = 'HIERARCH ESO QC BERV', # [km/s]
    fwhm     = 'HIERARCH ESO QC CCF FWHM', # CCF FWHM [km/s] 
    fwhm_err = 'HIERARCH ESO QC CCF FWHM ERROR', # Uncertainty on CCF FWHM [m/s?]
    cont     = 'HIERARCH ESO QC CCF CONTRAST', # CCF contrast [%]                
    cont_err = 'HIERARCH ESO QC CCF CONTRAST ERROR', # CCF contrast error [%] 
    bis      = 'HIERARCH ESO QC CCF BIS SPAN', # CCF bisector span [km/s]     
    bis_err  = 'HIERARCH ESO QC CCF BIS SPAN ERROR', # CCF bisector span err
)

# SP_HDRS = dict(
#     instr = 'INSTRUME',
#     obj = 'OBJECT',
#     bjd = 'HIERARCH ESO QC BJD',
#     date_obs = 'DATE-OBS',
#     exptime = 'EXPTIME',
#     snr60 = 'HIERARCH ESO QC ORDER85 SNR',
#     airmass_start = 'HIERARCH ESO TEL3 AIRM START',
#     airmass_end = 'HIERARCH ESO TEL3 AIRM END',
#     seeing_start = 'HIERARCH ESO TEL3 AMBI FWHM START', # [arcsec]
#     seeing_end = 'HIERARCH ESO TEL3 AMBI FWHM END', # [arcsec]
#     ra = 'RA',      # [deg]
#     dec = 'DEC',     # [deg]
#     prog_id = 'HIERARCH ESO OBS PROG ID',
#     pi_coi = 'PI-COI',
#     ftype = 'HIERARCH ESO PRO CATG',
#     drs_id = 'HIERARCH ESO PRO REC1 DRS ID',
#     targ_rv = 'HIERARCH ESO OCS OBJ RV',
#     rv_step = 'HIERARCH ESO RV STEP', #=     0.5 / Coordinate increment per pixel
#     rv_ref_pix = 'HIERARCH ESO RV START', #= -21.11 / Coordinate at reference pixel
#     berv = 'HIERARCH ESO QC BERV', # [km/s]
#     rv = 'HIERARCH ESO QC CCF RV', #= -1.19840358654132 / Radial velocity [km/s]        
#     rv_err = 'HIERARCH ESO QC CCF RV ERROR', #= 0.00053474023623152 / Uncertainty on radial veloc
#     fwhm = 'HIERARCH ESO QC CCF FWHM', #= 9.4897932752335 / CCF FWHM [km/s]                  
#     fwhm_err = 'HIERARCH ESO QC CCF FWHM ERROR', #= 0.00106948047246304 / Uncertainty on CCF FWHM [
#     cont = 'HIERARCH ESO QC CCF CONTRAST', #= 59.8806426078076 / CCF contrast %                
#     cont_err = 'HIERARCH ESO QC CCF CONTRAST ERROR', #= 0.00674842708267666 / CCF contrast error % 
#     continuum = 'HIERARCH ESO QC CCF CONTINUUM', #= 2560440.01971191 / CCF continuum level [e-]     
#     mask = 'HIERARCH ESO QC CCF MASK', #= 'F9      ' / CCF mask used                           
#     flux_asym = 'HIERARCH ESO QC CCF FLUX ASYMMETRY', #= 0.0126966418595723 / CCF asymmetry (km/s)  
#     flux_asym_err = 'HIERARCH ESO QC CCF FLUX ASYMMETRY ERROR', #= 0.000921109029304877 / CCF asymmetry 
#     bis = 'HIERARCH ESO QC CCF BIS SPAN', #= -0.01440933698031 / CCF bisector span (km/s)     
#     bis_err = 'HIERARCH ESO QC CCF BIS SPAN ERROR', #= 0.00106948047246304 / CCF bisector span err
#     blaze_file = "HIERARCH ESO PRO REC2 CAL12 NAME"
# )




class ESPRESSO:

    def __init__(self, file, hdu, obj_in=None, verb=False, save_spec=False, save_ccf=False, save_bis=False, output_path=None, add_spec_hdrs=None, add_ccf_hdrs=None, add_bis_hdrs=None, get_bis=True, get_ccf=True):
        """Reads ESPRESSO S1D and S2D fits files and extracts the spectrum and relevant header data.

        Args:
            file (str): Spectrum fits file with respective path.
            hdu ([type]): HDU of the spectrum fits file.
            obj_in (str): Identification of the target to override the target ID stored
                in the fits file.
            verb (bool, optional): If True prints information. Defaults to False.
            save_spec (bool, optional): If True saves the spectrum as a csv. Defaults to False.
            save_ccf (bool, optional): If True saves the CCF profile as a csv. Defaults to
                False.
            save_bis (bool, optional): If True saves the CCF bisector as a csv. Defaults to
                False.
            output_path (str, optional): Path for the output. Defaults to None.
            get_bis (bool, optional): Extract the CCF bisector and headers. Defaults to False.
            get_ccf (bool, optional): Extract the CCF profile and headers. Defaults to False.
            add_spec_hdrs (dict, optional): Dictionary with headers to be extracted from the
                main spectrum fits file. Defaults to None.
            add_ccf_hdrs (dict, optional): Dictionary with headers to be extracted from the
                'ccf' fits file. Defaults to None.
            add_bis_hdrs (dict, optional): Dictionary with headers to be extracted from the
                'bis' fits file. Defaults to None.
        """
        # Get additional headers:
        if add_spec_hdrs:
            spec_hdrs.update(add_spec_hdrs)
        
        flg = 'OK'
        instr = 'ESPRESSO'

        printif("Reading spectrum file", verb)

        # Create dictionary to hold the spectrum data
        spec = dict()

        # Get spectrum and file headers:
        spec, hdr = self._read_fits_espresso(hdu=hdu)

        # Get spectrum selected header values:
        headers = read_headers(hdr, spec_hdrs, data=None)

        # Get target:
        headers['obj'] = self._get_target(hdr, instr, verb)

        # Get median SNR:
        snr_med, _ = self._get_snr(hdr)
        headers['snr_med'] = snr_med

        # if no instrument keyword is present in fits file:
        if headers['instr'] is None:
            headers['instr'] = instr

        # Identify file type:
        if len(spec['flux_raw'].shape) == 1:
            headers['ftype'] = 'S1D'
        if len(spec['flux_raw'].shape) == 2:
            headers['ftype'] = 'S2D'

        for key in headers:
            if key in ['rv', 'fwhm', 'bis', 'berv']:
                headers[key] = headers[key]*1e3 # to m/s


        # shift wave to target rest frame:
        spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['rv'])

        spec['flux'] = spec['flux_raw']
        headers['spec_flg'] = 'OK'

        max_RON = 8 # [e-/pix] from the pipeline manual, pag. 16
        CONAD = 1.1 # [e-/ADU] from the pipeline manual, pag. 16
        headers['noise'] = max_RON * CONAD

        # CCF profile data:
        if get_ccf:
            printif("Reading CCF file", verb)
            ccf_file, _ = self._search_file(file, type='CCF', verb=False)

            hdu = fits.open(ccf_file)

            rv_step = hdr['HIERARCH ESO RV STEP'] # [km/s]
            rv_ref_pix = hdr['HIERARCH ESO RV START'] # [km/s]

            rv_grid = np.array([rv_ref_pix+ i * rv_step for i in range(hdu[1].data[0].size)])

            ccf_profile = dict(
                rv          = rv_grid,
                profile     = hdu[1].data,
                profile_err = hdu[2].data
            )

            self.ccf_profile = ccf_profile


        # save spectrum:
        if save_spec:
            if len(spec['wave'].shape) == 1:
                df = pd.DataFrame(spec)

            if len(spec['wave'].shape) == 2:
                data = {}
                for order in np.arange(spec['wave'].shape[0]):
                    for key in spec:
                        data[key + '_' + str(order)] = spec[key][order]
                df = pd.DataFrame(data)

            save_name = f"{headers['obj']}_{headers['instr']}_{headers['date_obs']}_spec.csv"
            self._save_object_data(df, headers['obj'], save_name, output_path, verb)

        # save CCF profile:
        if get_ccf and save_ccf:
            data = {}
            data['rv'] = ccf_profile['rv']
            for order in np.arange(ccf_profile['profile'].shape[0]):
                data['profile' + '_' + str(order)] = ccf_profile['profile'][order]
            df = pd.DataFrame(data)

            save_name = f"{headers['obj']}_{headers['instr']}_{headers['date_obs']}_ccf.csv"
            self._save_object_data(df, headers['obj'], save_name, output_path, verb)


        self.spectrum = spec
        self.headers = headers



    def _get_target(self, hdr, instr, verb):
        """
        Returns the object targeted in the fits file 'fits_file'.
        """
        if instr == 'HARPS':
            obs = 'ESO'
        if instr == 'HARPN':
            obs = 'TNG'

        try:
            obj = hdr['OBJECT']
        except KeyError:
            try:
                obj = hdr[f'{obs} OBS TARG NAME']
            except KeyError:
                printif("*** ERROR: Cannot identify object.", verb)
                return None

        return obj


    def _get_snr(self, hdr):
        i = 1 # espresso orders start at 1
        snr_orders = []
        while True:
            try:
                snr_orders.append(hdr[f'HIERARCH ESO QC ORDER{i} SNR'])
                i += 1
            except:
                break

        snr_orders = np.array(snr_orders)
        snr_med = np.median(snr_orders)
        
        return snr_med, snr_orders



    def _read_fits_espresso(self, file=None, hdu=None):

        if not hdu:
            hdu = fits.open(file)

        hdr = hdu[0].header

        if hdr['HIERARCH ESO PRO CATG'] == 'S1D_A':
            spec = dict(
                flux_raw = hdu[1].data['flux'], # debalzed
                flux_err = hdu[1].data['error'],
                wave_raw = hdu[1].data['wavelength_air'], # already in target rest frame
                #wave_raw = hdu[1].data['wavelength'], # already in target rest frame
            )

        elif hdr['HIERARCH ESO PRO CATG'] == 'S2D_A':
            spec = dict(
                flux_raw = hdu[1].data, # blazed
                flux_err = hdu[2].data,
                #wave_raw = hdu[4].data, # not in stellar rest frame
                wave_raw = hdu[5].data # wavelength air (bary)
            )
        
        else:
            raise ValueError(f"ERROR: File '{hdr['HIERARCH ESO PRO CATG']}' not implemented")

        hdu.close()

        return spec, hdr


    def _search_file(self, fits_file, type='ccf', verb=True):
        info = os.path.basename(fits_file).split('_')[0]
        path = os.path.dirname(fits_file)
        search_string = os.path.join(path, f"{info}_{type}_A.fits")

        file = glob.glob(search_string)

        flg = "OK"

        if len(file) == 0:
            err_msg = f"*** ERROR: {os.path.basename(__file__)}: No {type.upper()} file Found. Searched filename: {os.path.basename(search_string)}"
            flg = f"No{type.upper()}file"
            file = None

            printif(err_msg, verb)

            return file, flg

        elif len(file) == 1:
            file = file[0]

        elif len(file) > 1:
            err_msg = f"*** WARNING: Number of {type.upper()} files > 1. Using the first in the list. File used: {os.path.basename(file[0])}"

            printif(err_msg, verb)

            file = file[0]

        return file, flg


    def _save_object_data(self, df, obj, save_name, output_path, verb):
        if output_path:
            output = os.path.join(output_path, obj)
        else:
            output = obj

        if not os.path.isdir(output):
            os.makedirs(output)

        save_name = os.path.join(output, save_name)

        df.to_csv(save_name, index=False)
        printif(f"Data saved to {save_name}", verb)


