import sys, os
import numpy as np
import glob
import pandas as pd

from ._spec_tools import printif, read_fits, read_headers, wave_star_rest_frame, wave_corr_berv



spec_hdrs = dict(
    rv_lp       = 'HIERARCH ESO OBS TARG RADVEL', # low precision RV
    airmass     = 'HIERARCH ESO OBS TARG AIRMASS',
    obj         = 'HIERARCH ESO TEL TARG NAME',
    observer    = 'HIERARCH ESO OBS OBSERVER',
    date_obs    = 'HIERARCH ESO OBS FDATE',
    ra          = 'HIERARCH ESO TEL TARG ALPHA',
    dec         = 'HIERARCH ESO TEL TARG ALPHA',
    drs         = 'HIERARCH ESO DRS VERSION',
    instr       = 'HIERARCH ESO OBS INSTRUMENT',
    sigdet      = 'HIERARCH ESO DRS CCD SIGDET', # ccd readd out noise [e-]
    conad       = 'HIERARCH ESO DRS CCD CONAD', # CCD conv factor [e-/ADU]
    gain        = 'HIERARCH ESO CORA CCD GAIN', # CCD gain [e-/ADU]
    blaze_file  = 'HIERARCH ESO DRS BLAZE FILE',
    snr0        = 'HIERARCH ESO DRS SPE EXT SN0', # SNR order 0 (first)
    snr67       = 'HIERARCH ESO DRS SPE EXT SN67', # SNR order 67 (last)
    berv        = 'HIERARCH ESO DRS BERV', # Barycentric Earth Radial Velocity [km/s]
    wave_file   = 'HIERARCH ESO DRS CAL TH FILE',
    bjd         = 'HIERARCH ESO DRS BJD'
)


class CORALIE:
    def __init__(self, file, hdu, obj_in=None, verb=False, save_spec=False, out_path=None):

        flg = 'OK'
        instr = 'CORALIE'

        printif("Reading spectrum file", verb)

        # Create dictionary to hold the spectrum data
        spec = dict()

        # Get spectrum and file headers:
        spec['wave_raw'], spec['flux_raw'], hdr = read_fits(hdu=hdu, instr=instr)

        # Get spectrum selected header values:
        headers = read_headers(hdr, spec_hdrs, data=None)

        # Get target:
        #headers['obj'] = self._get_target(hdr, instr, verb)

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

        for key in headers.keys():
            if key in ['rv_lp', 'rv', 'berv']:
                headers[key] *= 1e3 # to m/s

        spec['wave'] = wave_star_rest_frame(spec['wave_raw'], headers['rv_lp'])

        # Flux photon noise:
        spec['flux_err'] = np.sqrt(abs(spec['flux_raw']))

        # Deblaze:
        if headers['ftype'] == 's1d':
            # Already deblazed
            spec['flux'] = spec['flux_raw']
            del spec['flux_raw']


        headers['noise'] = headers['sigdet'] * headers['conad'] # RON per pixel


        # output:
        self.spectrum = spec      # spectrum dict (must have 'wave' and 'flux')
        self.headers = headers    # all selected headers dict




    def _get_snr(self, hdr, instr):
        if instr in ['HARPS', 'CORALIE']:
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