# THIS FILE IS AN ARCHIVE TO BE DELETED


# Not used:
# HARPS
def get_headers():
    obs = 'ESO'

    spec_headers = dict(
        naxis    = 'NAXIS', # number of data axis
        naxis1   = 'NAXIS1', # length of data axis 1
        crpix1   = 'CRPIX1', # reference pixel
        cdelt1   = 'CDELT1', # Coordinate increment par pixel 
        crval1   = 'CRVAL1', # coordinate at reference pixel
        ctype1   = 'CTYPE1', # units of coordinate
        bunit    = 'BUNIT', # units of data values
        #tel      = 'TELESCOP',
        instr    = 'INSTRUME',
        date_obs = 'DATE-OBS',
        exptime  = 'EXPTIME',
        ra       = 'RA',
        dec      = 'DEC',
        bjd      = f'HIERARCH {obs} DRS BJD',
        snr7     = f'HIERARCH {obs} DRS SPE EXT SN{7}',
        snr50    = f'HIERARCH {obs} DRS SPE EXT SN{50}',
        sigdet   = f'HIERARCH {obs} DRS CCD SIGDET', #CCD Readout Noise [e-] 
        gain     = f'HIERARCH {obs} DRS CCD CONAD', #CCD conversion factor [e-/ADU]
        airmass_start = f'HIERARCH {obs} TEL AIRM START',
        airmass_end   = f'HIERARCH {obs} TEL AIRM END',
        seeing_start  = f'HIERARCH {obs} TEL AMBI FWHM START',
        seeing_end    = f'HIERARCH {obs} TEL AMBI FWHM END',
        catg     = f'HIERARCH {obs} DPR CATG',
        type     = f'HIERARCH {obs} DPR TYPE',
        blaze_file = f'HIERARCH {obs} DRS BLAZE FILE',
        prog_id  = f'HIERARCH {obs} OBS PROG ID',
        #pi_coi   = 'PI-COI',
        pi_coi   = f'HIERARCH {obs} OBS PI-COI NAME',
        drs      = f'HIERARCH {obs} DRS VERSION',
        targ_rv   = f'HIERARCH {obs} TEL TARG RADVEL',  # [km/s]
        #dome     = f'HIERARCH {obs} TEL DOME STATUS'
    )

    ccf_headers = dict(
        rv          = f"HIERARCH {obs} DRS CCF RVC",     # [km/s] (drift corrected)
        rv_err      = f"HIERARCH {obs} DRS DVRMS",       # [m/s]
        berv        = f"HIERARCH {obs} DRS BERV",        # [km/s]
        ccf_noise   = f'HIERARCH {obs} DRS CCF NOISE',   # [km/s] Photon noise on CCF RV
        fwhm        = f'HIERARCH {obs} DRS CCF FWHM',    # [km/s]
        cont        = f'HIERARCH {obs} DRS CCF CONTRAST', # [%]
        mask        = f"HIERARCH {obs} DRS CCF MASK"
    )

    bis_headers = dict(
        bis = f'HIERARCH {obs} DRS BIS SPAN' # [km/s]
    )

    return spec_headers, ccf_headers, bis_headers


# Not used:
def get_target(hdr):
    """
    Returns the object targeted in the fits file 'fits_file'.
    """
    try:
        obj = hdr['OBJECT']
    except KeyError:
        try:
            obj = hdr['ESO OBS TARG NAME']
        except KeyError:
            try:
                obj = hdr['TNG OBS TARG NAME']
            except KeyError:
                print("*** ERROR: Cannot identify object.")
                return None

    return obj


# SPIROU
spec_hdrs_0 = dict(
    obj     = 'OBJECT',
    instr   = 'INSTRUME',
    exptime = 'EXPTIME', # [sec]
    gain = hdr['GAIN'], # [e-/ADU]
    rdnoise = hdr['RDNOISE'], # [e-]
    pi_name = hdr['PI_NAME'],
    run_id = hdr['RUNID'],
    date_obs = hdr['DATE-OBS'] + 'T' + hdr['UTIME'],
    ra = hdr['RA_DEG'],
    dec = hdr['DEC_DEG'],
    airmass = hdr['AIRMASS'], # airmass at start
    drs_ver = hdr['VERSION'],
)