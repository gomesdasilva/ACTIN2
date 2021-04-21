#!/anaconda3/bin/python
"""
ACTIN 1.4: Line Parameters
"""
import sys, os
import numpy as np
import matplotlib.pylab as plt
import glob

# Local imports:
import fits_tools

import pandas as pd


# pandas version of get_config:
def get_ind_table(config_file=None):
    if config_file is None:
        file_name = "config_lines.txt"
        config_file = os.path.join(os.path.dirname(__file__), file_name)

    df = pd.read_csv(config_file, sep='\s+', comment='#', header=0)
    df = df.iloc[1:]

    return df


#get_ind_table()



def get_config(config_file=None):
    if config_file is None:
        config_file = "config_lines.txt"
    table = []
    with open(config_file) as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            if line.startswith('\n'):
                continue
            else:
                table.append(line)
    header = table[0].strip().split('\t')

    # remove empty keys
    header = [key for key in header if key]

    lines = []
    for line in table[2:]:
        lines.append(line.strip().split('\t'))

    ind_table = {}
    for key in header:
        ind_table[key] = []

    for line in lines:
        for i, (val, key) in enumerate(zip(line, header)):
            ind_table[key].append(line[i])

    indices = list(set(ind_table['ind_id']))

    return ind_table


def get_indices(ind_table):

    return list(set(ind_table['ind_id']))

def get_line_ids(index, ind_table):

    lines = []
    for i, ind in enumerate(ind_table['ind_id']):
        if ind == index:
            lines.append(ind_table['ln_id'][i])
    return lines


#* ind_table must be a dataframe:
def get_line_param(ind_table, param, line_id):

    if line_id not in ind_table.ln_id.values:
        print(f"*** ERROR: Line ID '{line_id}' not recognized")
        return np.nan

    value = ind_table[ind_table.ln_id == line_id][param].values

    #value = [ind_table[param][k] for k in range(len(ind_table[param])) if ind_table["ln_id"][k] == line_id]

    try:
        value = float(value[0])
    except: value = value[0]
    return value


# dependencies: get_line_param
def get_line_params(ind_table, line_id):
    #line_params = {}
    #for param in list(ind_table):
    #    line_params[param] = get_line_param(ind_table, param, line_id)
    #print(ind_table)

    return ind_table[ind_table["ln_id"] == line_id]


def get_wlims(ctr, win):
    wmin = ctr - win/2
    wmax = ctr + win/2
    return wmin, wmax


# dependencies: get_line_params
def plot_line(wave, flux, ind_table, line_id):

    line_par = get_line_params(ind_table, line_id)

    ctr = line_par['ln_ctr'].values[0]
    win = line_par['ln_win'].values[0]
    if line_par.bandtype.values[0] == 'sq':
        wmin = ctr - win/2
        wmax = ctr + win/2
    elif line_par.bandtype.values[0] == 'tri':
        wmin = ctr - win
        wmax = ctr + win

    plt.plot(wave, flux, 'k-', lw=0.7)
    plt.ylabel("Relative Flux")
    plt.xlabel("Wavelength [Ang]")
    plt.axvline(ctr, c='k')
    plt.axvline(wmin, c='k', ls='--')
    plt.axvline(wmax, c='k', ls='--')
    plt.xlim(wmin-16*win, wmax+16*win)
    plt.title(f"{line_id} line")
    plt.show()


# use with HARPS s1d, plots for all fits files:
def plot_line_star(star, line_id):
    config_file = "/Users/jgsilva/Astrophysics/Packages/_dev/actin14/config_lines.txt"
    #config_file = os.path.join(os.getcwd(), "config_lines.txt")
    stars_path = "/Volumes/GAMMA/Ambre/stars"
    s1d_files = glob.glob(os.path.join(stars_path, star, "*_s1d_A.fits"))

    for s1d_file in s1d_files[:]:
        wave, flux, flg, hdr = fits_tools.get_spec(s1d_file, instr='HARPS')
        #wave, flux, hdr, radvel, rv, flg_calib = fits_tools.calib_spec(s1d_file, star)
        #if flg_calib == 'NotCalib': continue
        #print(radvel, rv, flg_calib)

        config_table = get_config(config_file)

        plot_line(wave, flux, config_table, line_id)


#plot_line_star('*_del_Eri', 'Ha06')

# MAKE CALC_FLUX FUNCTION WITH FRAC PIXELS AND COMPARE TO THE INDEX WITHOUT FRAC!


# Ratio between core line and region outside bandpass:
def line_ratio(wave, flux, ind_table, line_id):
    line_par = get_line_params(ind_table, line_id)

    ctr = line_par['ln_ctr'].values[0]
    win = line_par['ln_win'].values[0]
    wmin = ctr - win/2
    wmax = ctr + win/2

    # average flux in continuum:
    cond_r = (wave > wmax) & (wave < wmax +win/2)
    cond_l = (wave > wmin-win/2) & (wave < wmin)
    flux_r = np.mean(flux[cond_r])
    flux_l = np.mean(flux[cond_l])
    wave_r = wave[cond_r][-1]-wave[cond_r][0]
    wave_l = wave[cond_l][-1]-wave[cond_l][0]

    flux_cont = np.mean([flux_r, flux_l])#/(wave_r+wave_l)
    #print(flux_cont)

    cond_line = (wave > wmin) & (wave < wmax)
    flux_line = np.mean(flux[cond_line])#/(wave[cond_line][-1]-wave[cond_line][0])
    #print(flux_line)

    print(flux_line/flux_cont)
    #sys.exit()
    plt.plot(wave, flux, 'k-', lw=0.7)
    plt.ylabel("Relative Flux")
    plt.xlabel("Wavelength [Ang]")
    plt.axvline(ctr, c='k')
    plt.axvline(wmin, c='k', ls='--')
    plt.axvline(wmax, c='k', ls='--')

    plt.axhline(flux_cont, c='k', ls='-')
    plt.axhline(flux_line, c='r', ls='-')

    plt.xlim(wmin-2*win, wmax+2*win)
    plt.title(f"{line_id} line")
    plt.show()

# IMPORTANT: not using fractional pixels!
def flux_phot_line(wave, flux, conf_table, line_id, sptype='FGK'):
    """Used to calculate Sphot for the calibration of R'HK. Calculates Fphot for spectral line 'line_id'. Fphot is the flux between the line core bandpass and the full line bandpass (the photospheric contribution to the line flux). Method used by Suárez Mascareño et al. (2015)."""

    line_par = get_line_params(conf_table, line_id)

    ctr = line_par['ln_ctr'].values[0]
    win = line_par['ln_win'].values[0]
    wmin = ctr - win/2
    wmax = ctr + win/2

    # line core bandpasses from Suarez Mascareno et al. (2015)
    if sptype == 'FGK':
        win_core = 0.7
    if sptype == 'M':
        win_core = 0.4

    wmin_core = ctr - win_core/2
    wmax_core = ctr + win_core/2

    cond_l = np.logical_and((wave > wmin), (wave < wmin_core))
    cond_r = np.logical_and((wave > wmax_core), (wave < wmax))
    cond = cond_l | cond_r

    flux_phot = np.sum(flux[cond])
    return flux_phot


# IMPORTANT: not using fractional pixels!
def flux_line(wave, flux, conf_table, line_id):

    line_par = get_line_params(conf_table, line_id)

    ctr = line_par['ln_ctr'].values[0]
    win = line_par['ln_win'].values[0]
    wmin = ctr - win/2
    wmax = ctr + win/2

    cond = np.logical_and((wave > wmin), (wave < wmax))

    F = np.sum(flux[cond])
    F_err = np.sqrt(np.sum(flux[cond]**2))
    return F, F_err

def flux_region(wave, flux, wmin, wmax):
    mask = (wave >= wmin) & (wave <= wmax)
    F = np.sum(flux[mask])
    F_err = np.sqrt(np.sum(flux[mask]**2))
    return F, F_err

def flux_mean_region(wave, flux, wmin, wmax):
    mask = (wave >= wmin) & (wave <= wmax)
    F = np.sum(flux[mask])/wave[mask].size
    F_err = np.sqrt(np.sum(flux[mask]**2))/wave[mask].size
    return F, F_err


# IMPORTANT: not using fractional pixels!
def flux_mean_line(wave, flux, conf_table, line_id):

    line_par = get_line_params(conf_table, line_id)

    ctr = line_par['ln_ctr'].values[0]
    win = line_par['ln_win'].values[0]
    wmin = ctr - win/2
    wmax = ctr + win/2

    cond = np.logical_and((wave >= wmin), (wave <= wmax))

    F_mean = np.mean(flux[cond])
    F_mean_err = np.sqrt(np.sum(flux[cond]**2))/wave[cond].size
    return F_mean, F_mean_err


# END
