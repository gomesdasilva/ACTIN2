import sys, os
import numpy as np
import matplotlib.pylab as plt
import glob
import tqdm
import time
import pandas as pd
from scipy.interpolate import interp1d
import astropy.io.fits as pyfits

# local:
import fits_tools
import line_params

# actin files:
import actin2_func as func

# temporary:
import decorators


#! IMPORTANT: Info related to flux errors see: Ramm+ (2021, pag. 6) [https://arxiv.org/pdf/2101.06844.pdf]

#! IMPORTANT INFO REGARDING READOUT NOISE TO BE ADDED TO ACTIN 2:
#! https://www.eso.org/~ohainaut/ccd/sn.html



#*-----------------------------------------------------------------
#* CONFIG:
#*-----------------------------------------------------------------
VERSION = "2"
CONFIG_FILE = "config_lines.txt"
OUTPUT_PATH = os.path.join(os.pardir, "output")

FILE_TYPES = ['s1d', 'e2ds', 'S1D', 'S2D']

INSTRUMENTS = ['HARPS', 'HARPN', 'ESPRESSO']

#* HARPS:
# HARPS_SP_ORDERS = 72
# HARPS_OBS       = 'ESO'

# HEADERS = dict(
#         tel      = 'TELESCOP',
#         instr    = 'INSTRUME',
#         date_obs = 'DATE-OBS',
#         bjd      = f'HIERARCH {HARPS_OBS} DRS BJD'
#     )
#*-----------------------------------------------------------------
IND_TABLE = line_params.get_ind_table()


def get_instrument(instr):
    if instr == 'HARPS':
        orders  = 72
        obs     = 'ESO'
    elif instr == 'HARPN':
        orders  = 69
        obs     = 'TNG'
    elif instr == 'ESPRESSO':
        orders  = 170
        obs     = 'ESO'
    else:
        raise Exception(f"Unrecognized {instr}. Available instruments are: {INSTRUMENTS}")

    headers = dict(
        tel      = 'TELESCOP',
        instr    = 'INSTRUME',
        date_obs = 'DATE-OBS',
        bjd      = f'HIERARCH {obs} DRS BJD',
        sigdet   = f'HIERARCH {obs} DRS CCD SIGDET', #CCD Readout Noise [e-] 
        gain     = f'HIERARCH {obs} DRS CCD CONAD', #CCD conversion factor [e-/ADU]
        snr6     = f'HIERARCH {obs} DRS SPE EXT SN{6}',
        snr50    = f'HIERARCH {obs} DRS SPE EXT SN{50}',
        airmass  = 'AIRMASS',
        exptime  = 'EXPTIME',
        ra       = 'RA',
        dec      = 'DEC',
        prog_id  = f'HIERARCH {obs} OBS PROG ID',
        #pi_coi   = 'PI-COI',
        pi_coi   = f'HIERARCH {obs} OBS PI-COI NAME',
        drs      = f'HIERARCH {obs} DRS VERSION',
        radvel   = f'HIERARCH {obs} TEL TARG RADVEL',  # [km/s]
        catg     = f'HIERARCH {obs} DPR CATG',
        type     = f'HIERARCH {obs} DPR TYPE'
        #dome     = f'HIERARCH {obs} TEL DOME STATUS'
    )

    ccf_headers = dict(
        rv        = f"HIERARCH {obs} DRS CCF RVC",     # [km/s]
        rv_err    = f"HIERARCH {obs} DRS DVRMS",       # [m/s]
        berv      = f"HIERARCH {obs} DRS BERV",        # [km/s]
        ccf_noise = f'HIERARCH {obs} DRS CCF NOISE',   # [km/s]
        fwhm      = f'HIERARCH {obs} DRS CCF FWHM',    # [km/s]
        cont      = f'HIERARCH {obs} DRS CCF CONTRAST' # [%]
    )

    bis_headers = dict(
        bis       = f'HIERARCH {obs} DRS BIS SPAN' # [km/s]
    )
    return orders, obs, headers, ccf_headers, bis_headers


def printif(string, verb=True, title=False):
    if verb:
        print(string)
        if title:
            print("-"*len(string))


def list_indices():
    return IND_TABLE.ind_id.drop_duplicates().values


def get_file_type(file, verb=True):
    """Identify file type. Identifiable files are stored in FILE_TYPES"""
    ftype = None
    for ft in FILE_TYPES:
        if ft in os.path.basename(file):
            ftype = ft
            break
    if not ftype:
        sys.exit("\n*** ERROR: Unidentified file type.")
    return ftype


def save(df, star, save_data, output_dir):
    if not save_data:
        return

    instr = df.instr.values[0]
    ftype = df.ftype.values[0]
    save_name = f"{star}_{instr}_{ftype}_actin.csv"

    if not output_dir:
        path = os.path.join(sys.path[0], "output", save_name)
    else:
        path = os.path.join(output_dir, save_name)

    if not os.path.isdir(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))

    df.to_csv(path, index=False)
    print(f"\nData saved to '{path}'")

#! NOT USED
def read_ccf_file(ccf_file, instr='HARPS'):
    if not ccf_file:
        rv = np.nan; rv_err = np.nan; berv = np.nan; ccf_noise = np.nan
        fwhm = np.nan; cont = np.nan
    else:
        hdr = fits_tools.get_hdr(ccf_file)
        rv, rv_err, berv = fits_tools.get_rv_hdr(hdr, instr=instr)
        ccf_noise = fits_tools.get_ccf_noise_hdr(hdr, instr=instr)
        fwhm = fits_tools.get_fwhm_hdr(hdr, instr=instr)
        cont = fits_tools.get_cont_hdr(hdr, instr=instr)
        rv *= 1000 # to m/s
        berv *= 1000 # to m/s
        ccf_noise *= 1000 # to m/s
        fwhm *= 1000 # to m/s

    ccf_data = dict(
        rv = rv,
        rv_err = rv_err,
        berv = berv,
        ccf_noise = ccf_noise,
        fwhm = fwhm,
        cont = cont
    )
    return ccf_data


def deblaze(file, flux, hdr, flg, OBS):
    try:
        blaze_file = hdr[f'HIERARCH {OBS} DRS BLAZE FILE']
        blaze_path = os.path.join(os.path.dirname(file), blaze_file)
        blaze = fits_tools.get_fits_data(blaze_path)
        flux = flux/blaze
    except:
        print("*** WARNING: No deblaze file found.")
        flg = "NoBlaze"
    return flux, flg


def override_obj(obj_name):
    """
    Override object name with name given in obj_name option.
    """
    if type(obj_name) is list and len(obj_name) == 1:
        return obj_name[0]
    elif type(obj_name) is list and len(obj_name) > 1:
        sys.exit("*** ERROR: obj_name requires only one name, more than one given.")
    else: return obj_name


def get_target(hdr):
    """
    Returns the object targeted in the fits file 'fits_file'.
    """
    try:
        obj = hdr['OBJECT']
    except:
        try:
            obj = hdr['ESO OBS TARG NAME']
        except:
            try:
                obj = hdr['TNG OBS TARG NAME']
            except:
                print("*** ERROR: Cannot identify object.")
                return None

    return obj


# read_fits:
def read_fits(fits_file, instr='HARPS'):
    if not fits_file:
        return None, None
    hdu = pyfits.open(fits_file)
    data = hdu[0].data
    hdr = hdu[0].header
    return data, hdr


def search_file(fits_file, type='ccf', verb=True):
    info = os.path.basename(fits_file).split('_')[0]
    path = os.path.dirname(fits_file)
    search_string = os.path.join(path, f"{info}_{type}_*_A.fits")

    type_fits = glob.glob(search_string)

    flg = "OK"

    if len(type_fits) == 0:
        err_msg = f"*** ERROR: {os.path.basename(__file__)}: No {type.upper()} file Found. Searched filename: {os.path.basename(search_string)}"
        flg = f"No{type.upper()}file"
        type_fits = None
        if verb:
            print(err_msg)
        return type_fits, flg

    elif len(type_fits) == 1:
        type_fits = type_fits[0]

    elif len(type_fits) > 1:
        err_msg = f"*** WARNING: Number of {type.upper()} files > 1. Using the first in the list. File used: {os.path.basename(type_fits[0])}"
        if verb:
            print(err_msg)
        type_fits = type_fits[0]

    return type_fits, flg

# Read fits header data usibg 'headers' dictionary:
def read_hdr_data(hdr, headers, data=None):
    """Read fits header data using 'headers' dictionary.
    Result is included in 'data' dictionary if not 'None'."""
    if not data:
        data = {}

    for key, hdr_id in zip(headers.keys(), headers.values()):
        if not hdr:
            data[key] = np.nan
            continue
        try:
            data[key] = hdr[hdr_id]
        except KeyError:
            data[key] = np.nan
    return data


def calib_wave_s1d(wave_raw, rv):
    """Correct wavelength doppler effect.
    rv must be in m/s"""
    c = 299792458.0 # [m/s]
    #rv *= 1000 # [m/s]

    # wavelength calibration:
    dwave = rv * wave_raw / c
    wave = wave_raw - dwave

    return wave


def read_spec(file, rv_in=None, instr='HARPS', verb=True):
    """
    flg : string
        OK: Everything is OK
        NoCCF: No CCF file found or no RV obtained from CCF
        RVin: Using external RV to calibrate wavelength
    """
    data = dict() # is returned as final output
    spec = dict() # not returned

    flg = 'OK'

    # TODO: Function to get 'instr' automaticaly
    orders, obs, headers, ccf_headers, bis_headers = get_instrument(instr)

    ftype = get_file_type(file, verb=verb)

    # FITS FILE:
    spec['flux'], hdr = read_fits(file, instr='HARPS')

    # Get object:
    data['obj'] = get_target(hdr)

    data['file']  = os.path.basename(file)
    data['ftype'] = ftype

    # Get data from fits headers:
    data = read_hdr_data(hdr, headers, data=data)

    snr_list = [hdr[f'HIERARCH {obs} DRS SPE EXT SN{k}'] for k in range(orders)]
    data['snr_med'] = np.median(snr_list)


    # CCF FILE:
    ccf_file, flg = search_file(file, type='ccf', verb=False)
    spec['ccf_profile'], ccf_hdr = read_fits(ccf_file, instr='HARPS')
    ccf_data = read_hdr_data(ccf_hdr, ccf_headers, data=None)

    # BIS FILE:
    bis_file, _ = search_file(file, type='bis', verb=False)
    spec['bisector'], bis_hdr = read_fits(bis_file, instr='HARPS')
    bis_data = read_hdr_data(bis_hdr, bis_headers, data=None)

    # Convert units:
    if not np.isnan(ccf_data['rv']):
        for key in ccf_data.keys():
            if key in ['rv', 'berv', 'ccf_noise', 'fwhm']:
                ccf_data[key] *= 1000 # convert to m/s
    else:
        flg = "NoCCF"
        # If no CCF RV use fits radvel:
        ccf_data['rv'] = data['radvel'] * 1000 # to m/s

    # Use external RV to calibrate wavelength
    if rv_in:
        flg = "RVin"
        ccf_data['rv'] = rv_in

    if ccf_file:
        data['ccf_file']  = os.path.basename(ccf_file)

    # Deblaze:
    if ftype == 's1d':
        spec['dflux'] = np.asarray([])
    if ftype == 'e2ds':
        spec['dflux'], flg = deblaze(file, spec['flux'], hdr, flg, obs)

    # Calculate and calibrate wavelength:
    if ftype == 's1d' and not np.isnan(ccf_data['rv']):
        spec['wave_raw'] = fits_tools.calc_wave_raw_s1d(hdr, instr='HARPS')
        spec['wave'] = calib_wave_s1d(spec['wave_raw'], ccf_data['rv'])
    elif ftype == 'e2ds' and not np.isnan(ccf_data['rv']):
        spec['wave_raw'] = fits_tools.calc_wave_raw_e2ds(hdr, spec['flux'].shape, obs)
        spec['wave'] = fits_tools.calib_wave_e2ds(spec['wave_raw'], ccf_data['rv'], ccf_data['berv'])
    else:
        spec = None

    #! This is not correct, should be n_pix instead of 109
    data['noise'] = data['sigdet'] * np.sqrt(109) * data['gain']
    #data['noise'] = 0.0

    data['flg']   = flg
    
    data.update(ccf_data)
    data.update(bis_data)

    # Round all float data in dictionary:
    def round_dict(data, digit):    
        for key, value in zip(data.keys(), data.values()):
            if isinstance(value, float):
                data[key] = round(value, digit)
            else:
                continue
        return data


    data = round_dict(data, 5)

    #print(data)

    return data, spec


def process_spec(file, index_id, rv_in=None, step=1e-5, verb=True):
    """Uses 'spec' from 'read_file' to calculate things then update 'data' dictionary. If no 'spec', returns original 'data'."""

    data, spec = read_spec(file, rv_in=rv_in, instr='HARPS', verb=verb)

    if not spec or not index_id:
        return data

    # TODO: Add section to compute co-added CCF profile

    # CALCULATE INDICES:
    flux = spec['flux']
    dflux = spec['dflux']
    wave = spec['wave']
    noise = data['noise']

    actin = {}
    for index in index_id:
        ind, ind_err, r_neg_flux = func.calc_index(wave, flux, dflux, noise, index, IND_TABLE, step=step, verb=verb)
   
        actin[index] = round(ind, 6)
        actin[index + "_err"] = round(ind_err, 6)
        actin[index + "_rneg"] = round(r_neg_flux, 4)

    if np.isnan(actin[index]):
        data['flg'] = "specERR"

    data.update(actin)


    return data


# WORKING:
@decorators.timeit
def run_folder(files, index_id=["I_CaII"], obj=None, rv_in=None, save_data=False, output_dir=None, verb=True, bar=True):
    printif("\n\nRUNNING ACTIN 2", verb=verb, title=True)

    # Get input arguments: - not used
    input_args = locals()

    # Check that indices are in IND_TABLE:
    if index_id:
        for index in index_id:
            if index not in IND_TABLE.ind_id.values:
                raise Exception(f"'{index}' not in 'config_file.txt'. Available indices are {list_indices()}.")

    # Enforce 'files' type as list:
    if not isinstance(files, (list, np.ndarray, tuple)):
        files = [files]

    files.sort()

    if bar:
        files = tqdm.tqdm(files)

    actin = [] 
    for file in files:
        time.sleep(0.01)
        data = process_spec(file, index_id, rv_in=rv_in, step=1e-5, verb=verb)
        if obj: data['obj'] = override_obj(obj)
        actin.append(data)


    # Create dataframe:
    df = pd.DataFrame(actin)

    # Save data:
    save(df, data['obj'], save_data, output_dir)

    return df
