import sys, os
import numpy as np
import glob
import pandas as pd

from ._spec_tools import printif, read_fits, read_headers, wave_star_rest_frame, wave_corr_berv

"""
IMPORTANT NOTE: RV error is quaadratic sum of:
HIERARCH ESO DRS CAL TH ERROR [m/s]
HIERARCH ESO DRS CCF NOISE [km/s]
HIERARCH ESO DRS DRIFT NOISE [m/s]
"""
obs = 'ESO'

spec_hdrs = dict(
    obj        = 'OBJECT',
    instr      = 'INSTRUME',
    date_obs   = 'DATE-OBS',
    bjd        = f'HIERARCH {obs} DRS BJD',
    drs        = f'HIERARCH {obs} DRS VERSION',
    exptime    = f'EXPTIME',
    ra       = 'RA',
    dec      = 'DEC',
    snr7     = f'HIERARCH {obs} DRS SPE EXT SN{7}',
    snr50    = f'HIERARCH {obs} DRS SPE EXT SN{50}',
    prog_id  = f'HIERARCH {obs} OBS PROG ID',
    pi_coi   = f'HIERARCH {obs} OBS PI-COI NAME',
    cal_th_err = f'HIERARCH {obs} DRS CAL TH ERROR', # estim err on wave sol(m/s)
)

ccf_hdrs = dict(
    rv          = f"HIERARCH {obs} DRS CCF RVC",     # [km/s] (drift corrected)
    dvrms       = f"HIERARCH {obs} DRS DVRMS",       # [m/s]
    berv        = f"HIERARCH {obs} DRS BERV",        # [km/s]
    ccf_noise   = f'HIERARCH {obs} DRS CCF NOISE',   # [km/s] Photon noise on CCF RV
    fwhm        = f'HIERARCH {obs} DRS CCF FWHM',    # [km/s]
    cont        = f'HIERARCH {obs} DRS CCF CONTRAST', # [%]
    mask        = f"HIERARCH {obs} DRS CCF MASK",
    drift_noise = f'HIERARCH {obs} DRS DRIFT NOISE' # Th Drift photon noise [m/s]
)

bis_hdrs = dict(
        bis = f'HIERARCH {obs} DRS BIS SPAN' # [km/s]
    )




class HARPS:

    def __init__(self, file, hdu, obj_in=None, verb=False, save_spec=False, save_ccf=False, save_bis=False, output_path=None, add_spec_hdrs=None, add_ccf_hdrs=None, add_bis_hdrs=None, get_bis=True, get_ccf=True):
        """Reads HARPS s1d and e2ds fits files and extracts the spectrum and relevant header data.

        Args:
            file ([type]): Spectrum fits file with respective path.
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
        if add_ccf_hdrs:
            ccf_hdrs.update(add_ccf_hdrs)
        if add_bis_hdrs:
            bis_hdrs.update(add_bis_hdrs)

        flg = 'OK'
        instr = 'HARPS'

        printif("Reading spectrum file", verb)

        # Create dictionary to hold the spectrum data
        spec = dict()

        # Get spectrum and file headers:
        spec['wave_raw'], spec['flux_raw'], hdr = read_fits(hdu=hdu, instr=instr)


        # Get spectrum selected header values:
        headers = read_headers(hdr, spec_hdrs, data=None)

        # Get target:
        headers['obj'] = self._get_target(hdr, instr, verb)

        # Get median SNR:
        snr_med, _ = self._get_snr(hdr, instr)
        headers['snr_med'] = snr_med

        # Add a different target name as 'obj_in':
        if obj_in:
            headers['obj_in'] = obj_in

        # if no instrument keyword is present in fits file:
        if headers['instr'] is None:
            headers['instr'] = instr

        # Identify file type:
        if len(spec['flux_raw'].shape) == 1:
            headers['ftype'] = 's1d'
        if len(spec['flux_raw'].shape) == 2:
            headers['ftype'] = 'e2ds'

        # CCF profile data:
        if get_ccf:
            printif("Reading CCF file", verb)
            ccf_file, _ = self._search_file(file, headers['ftype'], type='ccf', verb=verb)

        if get_ccf and ccf_file:
            ccf_profile = dict()

            ccf_profile['rv'], ccf_profile['profile'], ccf_hdr = read_fits(ccf_file, instr=instr)
            #ccf_profile['profile'] = ccf_profile['profile'][-1] # co-added profile

            headers = read_headers(ccf_hdr, ccf_hdrs, data=headers)

            for key in headers.keys():
                if key in ['rv', 'berv', 'ccf_noise', 'fwhm']:
                    headers[key] *= 1e3 # to m/s

            self.ccf_profile = ccf_profile

        # CCF bisector data:
        if get_bis:
            printif("Reading BIS file", verb)
            bis_file, _ = self._search_file(file, headers['ftype'], type='bis', verb=verb)

        if get_bis and bis_file:
            ccf_bisector = dict()

            ccf_bisector['bisector'], ccf_bisector['rv'], bis_hdr = read_fits(bis_file, instr=instr)

            headers = read_headers(bis_hdr, bis_hdrs, data=headers)

            self.ccf_bisector = ccf_bisector


        # Calibrate wavelength to stellar rest frame:
        if get_ccf and ccf_file:
            rv = headers['rv']
            if headers['ftype'] == 's1d' and not np.isnan(rv):
                spec['wave'] = wave_star_rest_frame(spec['wave_raw'], rv)

            elif headers['ftype'] == 'e2ds' and not np.isnan(rv) and 'berv' in headers:
                spec['wave'] = wave_corr_berv(spec['wave_raw'], headers['berv'])
                spec['wave'] = wave_star_rest_frame(spec['wave'], rv)
        else:
            rv = np.nan
            flg = "WaveNotCorr"
            printif("*** WARNING: Cannot shift spectrum to target rest frame", verb)

        # Flux photon noise:
        spec['flux_err'] = np.sqrt(abs(spec['flux_raw']))

        # Deblaze:
        if headers['ftype'] == 's1d':
            # Already deblazed
            spec['flux'] = spec['flux_raw']
            del spec['flux_raw']

        if headers['ftype'] == 'e2ds':
            blaze_file = hdr['HIERARCH ESO DRS BLAZE FILE']
            spec['flux'], flg = self._deblaze(file, spec['flux_raw'], blaze_file, flg, instr)


        sigdet = hdr['HIERARCH ESO DRS CCD SIGDET'] #CCD Readout Noise [e-] 
        gain = hdr['HIERARCH ESO DRS CCD CONAD']  #CCD conversion factor [e-/ADU]
        headers['noise'] = sigdet * gain # Read-Out-Noise per pixel

        headers['rv_err'] = np.sqrt(headers['ccf_noise']**2 + headers['drift_noise']**2 + headers['cal_th_err']**2)

        headers['spec_flg'] = flg


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

            #df = pd.DataFrame(ccf_profile)
            save_name = f"{headers['obj']}_{headers['instr']}_{headers['date_obs']}_ccf.csv"
            self._save_object_data(df, headers['obj'], save_name, output_path, verb)

        # save CCF bisector:
        if get_bis and save_bis:
            df = pd.DataFrame(ccf_bisector)
            save_name = f"{headers['obj']}_{headers['instr']}_{headers['date_obs']}_bis.csv"
            self._save_object_data(df, headers['obj'], save_name, output_path, verb)


        # output:
        self.spectrum = spec      # spectrum dict (must have 'wave' and 'flux')
        self.headers = headers    # all selected headers dict


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


    def _get_snr(self, hdr, instr):
        if instr == 'HARPS':
            obs = 'ESO'
        if instr == 'HARPN':
            obs = 'TNG'
        
        i = 0
        snr_orders = []
        while True:
            try:
                snr_orders.append(hdr[f'HIERARCH {obs} DRS SPE EXT SN{i}'])
                i += 1
            except:
                break

        snr_orders = np.array(snr_orders)
        snr_med = np.median(snr_orders)
        
        return snr_med, snr_orders


    def _search_file(self, fits_file, ftype, type='ccf', verb=True):
        info = os.path.basename(fits_file).split(f'_{ftype}')[0]
        path = os.path.dirname(fits_file)
        search_string = os.path.join(path, f"{info}_{type}_*_A.fits")

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


    def _deblaze(self, file, flux, blaze_file, flg, instr):
        try:
            blaze_path = os.path.join(os.path.dirname(file), blaze_file)
            blaze, _ = read_fits(blaze_path, instr=instr, calc_x=False)
            flux = flux/blaze
        except:
            flg = "NoBlaze"
        return flux, flg


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


