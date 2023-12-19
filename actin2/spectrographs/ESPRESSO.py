import sys, os
import numpy as np
import glob
import pandas as pd
from astropy.io import fits

from ._spec_tools import printif, read_headers, wave_star_rest_frame



spec_hdrs = dict(
    obj      = 'OBJECT',
    instr    = 'INSTRUME',
    date_obs = 'DATE-OBS',
    bjd      = 'HIERARCH ESO QC BJD',
    rv       = 'HIERARCH ESO QC CCF RV', # Radial velocity [km/s]
    rv_err   = 'HIERARCH ESO QC CCF RV ERROR', # Uncertainty on radial velocity [km/s]
    berv     = 'HIERARCH ESO QC BERV', # [km/s]
    fwhm     = 'HIERARCH ESO QC CCF FWHM', # CCF FWHM [km/s] 
    fwhm_err = 'HIERARCH ESO QC CCF FWHM ERROR', # Uncertainty on CCF FWHM [km/s]
    cont     = 'HIERARCH ESO QC CCF CONTRAST', # CCF contrast [%]                
    cont_err = 'HIERARCH ESO QC CCF CONTRAST ERROR', # CCF contrast error [%] 
    bis      = 'HIERARCH ESO QC CCF BIS SPAN', # CCF bisector span [km/s]     
    bis_err  = 'HIERARCH ESO QC CCF BIS SPAN ERROR', # CCF bisector span err
    exptime  = 'EXPTIME',
    snr50    = 'HIERARCH ESO QC ORDER50 SNR',
    spec_rv  = 'HIERARCH ESO OCS OBJ RV',
    ra       = 'RA',      # [deg]
    dec      = 'DEC',     # [deg]
    prog_id  = 'HIERARCH ESO OBS PROG ID',
    pi_coi   = 'PI-COI',
    ftype    = 'HIERARCH ESO PRO CATG',
    drs_id   = 'HIERARCH ESO PRO REC1 DRS ID',
    ccf_mask = 'HIERARCH ESO QC CCF MASK',

)

ccf_hdrs = dict()

bis_hdrs = dict()



class ESPRESSO:

    def __init__(self, hdu, file=None, rv_in="ccf", berv_in=None, verb=False, save_spec=False, save_ccf=False, save_bis=False, output_path=None, add_spec_hdrs=None, add_ccf_hdrs=None, add_bis_hdrs=None, get_bis=True, get_ccf=True):
        """Reads ESPRESSO S1D and S2D fits files and extracts the spectrum and relevant header data.

        Args:
            file (str): Spectrum fits file with respective path.
            hdu (fits HDU): HDU of the spectrum fits file.
            rv_in (str, float): RV to calibrate wavelength to stellar rest frame. If 'ccf' (default) use the CCF RV, if 'spec' use the spectrum headers RV, if float use input RV (in km/s).
            berv_in (None, float): Input BERV value to correct wavelength (m/s). If 'None' use HARPS DRS BERV.
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
        # Get additional headers for each file type:
        if add_spec_hdrs and isinstance(add_spec_hdrs, dict):
            spec_hdrs.update(add_spec_hdrs)
        if add_ccf_hdrs and isinstance(add_ccf_hdrs, dict):
            ccf_hdrs.update(add_ccf_hdrs)
        if add_bis_hdrs and isinstance(add_bis_hdrs, dict):
            bis_hdrs.update(add_bis_hdrs)
        
        flg = 'OK'
        instr = 'ESPRESSO'

        np.seterr(divide='ignore', invalid='ignore')

        self.file = file

        printif(f"Reading spectrum file: {file}", verb)

        # Create dictionary to hold the spectrum data
        spec = dict()

        # Get spectrum and file headers:
        spec, hdr = self._read_fits_espresso(hdu=hdu)

        # Get spectrum selected header values:
        headers = read_headers(hdr, spec_hdrs, data=None, verb=verb)

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
            
            spec['flux'] = spec['flux_raw']/spec['dlldata']
            flux_rel_err = spec['flux_err']/abs(spec['flux_raw'])
            spec['flux_err'] = flux_rel_err * abs(spec['flux'])


        # The RV used to calibrate the wavelength to the stellar rest frame can come from the CCF, the spectrum headers or as input:
        if rv_in == 'ccf':
            headers['rv_wave_corr'] = headers['rv']
            headers['rv_flg'] = 'CCF'

        elif rv_in == 'spec' and 'spec_rv' in headers:
            headers['rv_wave_corr'] = headers['spec_rv']
            headers['rv_flg'] = 'SPEC'

        elif isinstance(rv_in, float):
            headers['rv_wave_corr'] = rv_in # km/s
            headers['rv_flg'] = 'INPUT'

        else:
            print(f"{__class__.__name__} WARNING: could not calibrate wavelength to stellar rest frame, 'rv_in' input was '{rv_in}'")
            headers['rv_flg'] = 'NoRV'
            headers['rv_wave_corr'] = None
            self.spectrum = spec
            self.headers = headers
            return


        keys = ['rv', 'rv_err', 'berv', 'ccf_noise', 'fwhm', 'fwhm_err', 'bis', 'bis_err', 'spec_rv', 'rv_wave_corr']
        for key in keys:
            if key in headers:
                try:
                    headers[key] *= 1e3 # to m/s
                except TypeError:
                    continue


        # shift wave to target rest frame:
        spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['rv_wave_corr'])

        # flux already deblazed in S1D files
        if headers['ftype'] == 'S1D':
            spec['flux'] = spec['flux_raw']

        headers['spec_flg'] = 'OK'


        # CCF profile data:
        if get_ccf:
            printif("Reading CCF file", verb)
            ccf_file, _ = self._search_file(file, type='CCF', verb=False)

            if ccf_file is not None:
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
                printif(f"{__class__.__name__} ERROR: Cannot identify object.", verb)
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
                wave_raw = hdu[1].data['wavelength_air'],
            )

        elif hdr['HIERARCH ESO PRO CATG'] == 'S2D_A':
            spec = dict(
                flux_raw = hdu[1].data, # deblazed
                flux_err = hdu[2].data,
                wave_raw = hdu[5].data, # wavelength air (bary)
                dlldata = hdu[7].data, # for normalization
            )

            if file is not None:
                blaze_file = file.split("_S2D")[0] + "_S2D_BLAZE_A.fits"

                if not os.path.isfile(blaze_file):
                    printif(f"Blaze file {blaze_file} not found. Spectrum was not deblazed.")
                    return spec, hdr

                hdu = fits.open(blaze_file)
                spec['blaze'] = hdu[1].data
            else:
                blaze_file = "r." + hdu[0].header["ARCFILE"].split(".fits")[0] + "_S2D_BLAZE_A.fits"
                blaze_file = os.path.join(os.path.dirname(self.file), blaze_file)

                if not os.path.isfile(blaze_file):
                    printif(f"Blaze file {blaze_file} not found.")
                    return spec, hdr

                hdu = fits.open(blaze_file)
                spec['blaze'] = hdu[1].data
        
        else:
            raise ValueError(f"{__class__.__name__} ERROR: File '{hdr['HIERARCH ESO PRO CATG']}' not implemented")

        hdu.close()

        return spec, hdr


    def _search_file(self, fits_file, type='ccf', verb=True):
        info = os.path.basename(fits_file).split('_')[0]
        path = os.path.dirname(fits_file)
        search_string = os.path.join(path, f"{info}_{type}_A.fits")

        file = glob.glob(search_string)

        flg = "OK"

        if len(file) == 0:
            err_msg = f"{__class__.__name__} ERROR: {os.path.basename(__file__)}: No {type.upper()} file Found. Searched filename: {os.path.basename(search_string)}"
            flg = f"No{type.upper()}file"
            file = None

            printif(err_msg, verb)

            return file, flg

        elif len(file) == 1:
            file = file[0]

        elif len(file) > 1:
            err_msg = f"{__class__.__name__} WARNING: Number of {type.upper()} files > 1. Using the first in the list. File used: {os.path.basename(file[0])}"

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


