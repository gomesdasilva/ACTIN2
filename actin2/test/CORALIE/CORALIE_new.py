import numpy as np


spec_hdrs = dict(
    instr   = 'HIERARCH ESO OBS INSTRUMENT',
    bjd     = 'HIERARCH ESO DRS BJD',
    spec_rv = 'HIERARCH ESO OBS TARG RADVEL', # low precision RV [km/s]
    berv    = 'HIERARCH ESO DRS BERV', # Barycentric Earth Radial Velocity [km/s]    
)


class CORALIE_new:
    def __init__(self, hdu, **spec_kw):

        instr = 'CORALIE'

        # Create dictionary to hold the spectrum data
        spec = dict()

        # Obtain the flux and headers from fits HDU
        spec['flux_raw'] = hdu[0].data
        hdr = hdu[0].header
        hdu.close()

        # Calculate wavelength grid
        spec['wave_raw'] = hdr['CRVAL1'] + hdr['CDELT1'] * np.arange(hdr['NAXIS1'])

        # Get spectrum selected header values:        
        headers = {}
        for key, hdr_id in zip(spec_hdrs.keys(), spec_hdrs.values()):
            try:
                headers[key] = hdr[hdr_id]
            except KeyError:
                headers[key] = None

        headers['instr'] = instr

        # Convert RV and BERV to m/s
        for key in headers.keys():
            if key in ['spec_rv', 'berv']:
                headers[key] *= 1e3 # to m/s

        # Correct spectrum to stellar rest frame
        c = 299792458.0 # light velocity [m/s]
        dwave = headers['spec_rv'] * spec['wave_raw'] / c
        spec['wave'] = spec['wave_raw'] - dwave

        # Flux photon noise
        spec['flux_err'] = np.sqrt(abs(spec['flux_raw']))

        # s1d files already deblazed
        spec['flux'] = spec['flux_raw']


        # output:
        self.spectrum = spec      # spectrum dict (must have 'wave' and 'flux')
        self.headers = headers    # all selected headers dict


