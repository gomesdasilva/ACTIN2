
import sys, os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from scipy.interpolate import interp1d

# local:
import line_params as lp


def get_line_ids(index_id, ind_table):
    return ind_table[ind_table.ind_id == index_id].ln_id.values



def check_negflux(flux, verb=True):
    flux = np.asarray(flux)

    # Total absolute flux in given spectral line:
    tot_flux = np.sum(abs(flux))

    # Number of pixels with negative flux
    neg_pixs = flux[flux < 0].size

    # Negative flux in given spectral line
    neg_flux = np.sum(flux[flux < 0])

    # Positive flux in given spectral line
    pos_flux = np.sum(flux[flux > 0])

    # Negative flux ratio for a given spectral line:
    r_neg_ln = abs(neg_flux)/tot_flux

    flg_negflux = "OK"

    if neg_pixs > 0:
        if verb:
            print(f"Negative flux detected")
        flg_negflux = "negFlux"

    return r_neg_ln, neg_flux, tot_flux, flg_negflux



#* Interpolating between limiting pixels only
def inter_flux_band(wave, flux, dflux, noise, ctr, win, bptype='tri', step=1e-4, make_plot=False):
    """Calculate average flux in bandpass"""

    #res_ratio = np.average(np.diff(wave))/step

    if bptype == "tri":
        wmin = ctr - win; wmax = ctr + win
    if bptype == "sq":
        wmin = ctr - win/2; wmax = ctr + win/2

    def interpolate_band_lims(array, wave, wmin, wmax, step):
        res_ratio = np.average(np.diff(wave))/step

        mask = (wave >= wmin) & (wave <= wmax)
        wave_int_low = (wave[wave < wmin][-1], wave[wave >= wmin][0])
        wave_int_high = (wave[wave <= wmax][-1], wave[wave > wmax][0])

        interv_low = (wave >= wave_int_low[0]) & (wave <= wave_int_low[1])
        interv_high = (wave >= wave_int_high[0]) & (wave <= wave_int_high[1])

        wave_low = wave[interv_low]
        wave_high = wave[interv_high]

        array_low = array[interv_low]
        array_high = array[interv_high]

        interp_low = interp1d(wave_low, array_low, kind='linear')
        interp_high= interp1d(wave_high, array_high, kind='linear')

        wave_i_low = np.arange(min(wave_low), max(wave_low), step)
        array_i_low = interp_low(wave_i_low)

        wave_i_high = np.arange(min(wave_high), max(wave_high), step)
        array_i_high = interp_high(wave_i_high)

        wave_i = np.r_[wave_i_low, wave[mask], wave_i_high]

        array_i = np.r_[array_i_low/(array_i_low.size*res_ratio), array[mask], array_i_high/(array_i_high.size*res_ratio)]

        return wave_i, array_i


    wave_i, flux_i = interpolate_band_lims(flux, wave, wmin, wmax, step)

    if dflux.size != 0:
        _, dflux_i = interpolate_band_lims(dflux, wave, wmin, wmax, step)
    else:
        dflux_i = flux_i

    r_neg_ln, _, _, _ = check_negflux(dflux_i, verb=False)
    
    if bptype == 'tri':
        bp_i = 1 - np.abs(wave_i - ctr)/win
        bp_i = np.where(bp_i > 0, bp_i, bp_i*0.0)
    elif bptype == 'sq':
        bp_mask = (wave_i >= ctr - win/2) & (wave_i <= ctr + win/2)
        bp_i = np.where(bp_mask, 1, 0.0)

    if make_plot:
        plt.plot(wave_i, dflux_i, 'k.-', lw=0.7)
        plt.plot(wave_i, dflux_i.max()*bp_i, 'b-', lw=0.7)
        plt.axvline(wmin, c='r')
        plt.axvline(wmax, c='r')
        plt.xlabel(r"$\lambda$ [$\AA$]")
        plt.ylabel("Norm. Flux")
        plt.show()

    flux_band = sum(dflux_i * bp_i)/win
    flux_band_var= sum((flux_i + noise**2) * bp_i**2)/win**2

    return flux_band, np.sqrt(flux_band_var), r_neg_ln





def spec_order(wave_2d, flux_2d, dflux_2d, ln_ctr, ln_win, bandtype):
    """Choose best spectral order for line and returns wave and flux at order.
    """
    if bandtype == 'tri':
        wmin = ln_ctr - ln_win
        wmax = ln_ctr + ln_win
    if bandtype == 'sq':
        wmin = ln_ctr - ln_win/2
        wmax = ln_ctr + ln_win/2

    orders = []
    min_dist = []
    for i, wave in enumerate(wave_2d):
        if wmin > wave[1] and wmax < wave[-1]:
            orders.append(i)
            dist = ((wmin - wave_2d[i][0], wave_2d[i][-1] - wmax))
            min_dist.append(np.min(dist))

    # Selects the order with the higher distance to the wave edges:
    order = orders[np.argmax(min_dist)]

    if dflux_2d.size == 0:
        return wave_2d[order], flux_2d[order], dflux_2d

    return wave_2d[order], flux_2d[order], dflux_2d[order]



def calc_index(wave, flux, dflux, noise, index_id, ind_table, step=1e-4, verb=True):
    """Calculate activity index 'index_id' using lines data from 'ind_table'."""
    line_ids = get_line_ids(index_id, ind_table)

    # function to calculate flux:
    get_flux = inter_flux_band

    if not noise:
        noise = 0.0

    F_num = []
    F_num_err = []
    F_denum = []
    F_denum_err = []
    r_neg_flux = []
    for line in line_ids:
        #print("LINE:", line)
        df = lp.get_line_params(ind_table, line)

        ctr = df.ln_ctr.values[0]
        win = df.ln_win.values[0]
        bandtype = df.bandtype.values[0]
        ind_var = df.ind_var.values[0]

        if len(wave.shape) == 2:
            w, f, debf = spec_order(wave, flux, dflux, ctr, win, bandtype)
        if len(wave.shape) == 1:
            w = wave; f = flux; debf = dflux

        try:
            f, err, r_neg_ln = get_flux(w, f, debf, noise, ctr, win, bptype=bandtype, step=step, make_plot=False)
        except:
            return np.nan, np.nan, np.nan

        r_neg_flux.append(r_neg_ln)

        if ind_var.startswith("L"):
            F_num.append(f)
            F_num_err.append(err)

        elif ind_var.startswith("R"):
            F_denum.append(f)
            F_denum_err.append(err)

    F_num = np.asarray(F_num)
    F_denum = np.asarray(F_denum)
    F_num_err = np.asarray(F_num_err)
    F_denum_err = np.asarray(F_denum_err)

    if not F_denum.all():
        F_denum = 1.

    index = np.sum(F_num)/np.sum(F_denum)
    index_err = np.sqrt(np.sum(F_num_err**2) + index**2 * np.sum(F_denum_err**2)) / (np.sum(F_denum))
    r_neg_flux = max(r_neg_flux)#/len(r_neg_flux)
    return index, index_err, r_neg_flux




