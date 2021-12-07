import sys, os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from scipy.interpolate import interp1d


from .indtable import IndTable
from .spectrographs._spec_tools import printif


def spec_order(wave_2d, flux_2d, flux_2d_err, ln_ctr, ln_win, bandtype, type='dist', show_orders=False, verb=False):
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
        if wmin > wave[0] and wmax < wave[-1]:
            orders.append(i)
            dist = ((wmin - wave_2d[i][0], wave_2d[i][-1] - wmax))
            min_dist.append(np.min(dist))


    printif(f"available orders: {orders}", verb)

    if show_orders:
        plt.title("spec_order: orders available for selected line")
        for order in orders:
            plt.plot(wave_2d[order], flux_2d[order], ls='-', marker='', label=order)
        plt.axvline(ln_ctr, c='k', ls='-')
        plt.axvline(wmin, c='k', ls=':')
        plt.axvline(wmax, c='k', ls=':')
        plt.legend()
        plt.show()

    # Selects the order with the higher distance to the wave edges from the bandpass limits:
    if type == 'dist':
        order = orders[np.argmax(min_dist)]
        printif(f"selected: {order}", verb)

    return wave_2d[order], flux_2d[order], flux_2d_err[order]


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




class CalcIndices:

    # TODO: add the above functions to the class?
    # TODO: Try to separate total noise from photon noise and RON

    def __init__(self, spec, indices, step=1e-6, table_df=None, plot_lines=False, full_output=False, interp=True, verb=False):

        data = spec.headers

        try:
            flux = spec.spectrum['flux']
        except KeyError:
            printif("*** ERROR: There is no 'flux' keyword in the spectrum dictionary", verb)
            for index in indices:
                data[index] = np.nan
                data[index + "_err"] = np.nan
            self.data = data
        else:
            try:
                wave = spec.spectrum['wave']
            except KeyError:
                printif("*** ERROR: There is no 'wave' keyword in the spectrum dictionary", verb)
                for index in indices:
                    data[index] = np.nan
                    data[index + "_err"] = np.nan
                self.data = data

            else:

                try:
                    flux_err = spec.spectrum['flux_err']
                except KeyError:
                    flux_err = np.sqrt(abs(flux))

                try:
                    noise = spec.headers['noise']
                except KeyError:
                    noise = 0.0


                for index in indices:
                    printif(f" Computing {index}", verb=verb)
                    index_dict = self.calc_index(wave, flux, flux_err, noise, index, step, table_df, plot_lines, full_output, interp, verb)

                    data.update(index_dict)

                self.data = data


    def calc_index(self, wave, flux, flux_err, noise, index, step, table_df, plot_lines, full_output, interp, verb):

        printif("Running CalcIndex", verb)

        if table_df is None:
            ind_table = IndTable().table
        else:
            ind_table = table_df

        assert index in ind_table.ind_id.values, f"*** ERROR: {index} is not available in indices table"

        # get index from ind_table:
        df = ind_table[ind_table.ind_id == index]
        ln_ids = df.ln_id.values


        F_num = []
        F_num_err = []
        F_denum = []
        F_denum_err = []
        r_neg_flux = []
        lines = {}
        for line in ln_ids:
            mask = (df.ln_id == line)

            ctr = df.ln_ctr[mask].values[0]
            win = df.ln_win[mask].values[0]
            bandtype = df.bandtype[mask].values[0]
            const = df.ln_c[mask].values[0]
            ind_var = df.ind_var[mask].values[0]


            if len(wave.shape) == 2:
                w, f, f_err = spec_order(wave, flux, flux_err, ctr, win, bandtype)
            if len(wave.shape) == 1:
                w = wave; f = flux; f_err= flux_err

            
            F, F_err, rneg_ln = self.inter_flux_band(w, f, f_err, noise, ctr, win, bandtype, step, plot_lines, interp)

            lines[line + '_F'] = F
            lines[line + '_F_err'] = F_err
            lines[line + 'rneg'] = rneg_ln
                

            r_neg_flux.append(rneg_ln)

            if ind_var.startswith("L"):
                F_num.append(const * F)
                F_num_err.append(F_err)

            elif ind_var.startswith("R"):
                F_denum.append(const * F)
                F_denum_err.append(F_err)


        F_num = np.asarray(F_num)
        F_denum = np.asarray(F_denum)
        F_num_err = np.asarray(F_num_err)
        F_denum_err = np.asarray(F_denum_err)

        if not F_denum.all():
            F_denum = 1.

        val = np.sum(F_num)/np.sum(F_denum)

        val_err = np.sqrt(np.sum(F_num_err**2) + val**2 * np.sum(F_denum_err**2)) / (np.sum(F_denum))

        if np.all(r_neg_flux) == np.zeros_like(r_neg_flux).all():
            r_neg_flux = 0.0
        else:
            r_neg_flux = np.mean(r_neg_flux) # mean ratio of negative flux to total flux in bands


        data = {}
        data[index] = val
        data[index + '_err'] = val_err
        data[index + '_mrneg'] = r_neg_flux

        if full_output:
            data.update(lines)

        return data


    
    def inter_flux_band(self, wave, flux, flux_err, noise, ctr, win, bandtype, step, show_plot=False, interp=True):
        """Calculate average flux in bandpass.
        
        Interpolates between the limiting pixels and the bandpass limits to deal with the finite spectrograph resolution."""

        if bandtype == "tri":
            wmin = ctr - win; wmax = ctr + win
        if bandtype == "sq":
            wmin = ctr - win/2; wmax = ctr + win/2

        def interpolate_band_lims(array, wave, wmin, wmax, step):
            """Interpolate the flux between the limiting pixels and the bandpass limits wmin and wmax to reduce the effect of the finite spectrograph resolution.

            Args:
                array (ndarray): y-axis data.
                wave (ndarray): x-axis data.
                wmin (ndarray): Left limit of bandpass.
                wmax (ndarray): Right limit of bandpass.
                step (ndarray): Step used in the interpolation.

            Returns:
                ndarray: Interpolated x-axis data.
                ndarray: Interpolated y-axis data.
            """
            # res_ratio = np.average(np.diff(wave))/step

            mask = (wave >= wmin) & (wave <= wmax)

            # wavelength right before and after the window limits:
            wave_int_low = (wave[wave < wmin][-1], wave[wave >= wmin][0])
            wave_int_high = (wave[wave <= wmax][-1], wave[wave > wmax][0])

            interv_low = (wave >= wave_int_low[0]) & (wave <= wave_int_low[1])
            interv_high = (wave >= wave_int_high[0]) & (wave <= wave_int_high[1])

            wave_low = wave[interv_low]
            wave_high = wave[interv_high]

            array_low = array[interv_low]
            array_high = array[interv_high]

            interp_low = interp1d(wave_low, array_low, kind='linear', fill_value="extrapolate")
            interp_high= interp1d(wave_high, array_high, kind='linear', fill_value="extrapolate")


            wave_i_low = np.arange(min(wave_low), max(wave_low) + step, step)
            array_i_low = interp_low(wave_i_low)

            wave_i_high = np.arange(min(wave_high), max(wave_high) + step, step)
            array_i_high = interp_high(wave_i_high)

            # wave_i = np.r_[wave_i_low, wave[mask], wave_i_high]
            # array_i = np.r_[array_i_low/res_ratio, array[mask], array_i_high/res_ratio]

            # mask = (wave_i >= wmin) & (wave_i <= wmax)
            # array_i = array_i[mask]
            # wave_i = wave_i[mask]

            mask_low = (wave_i_low >= wmin)
            mask_high = (wave_i_high <= wmax)

            wave_i = np.r_[min(wave_i_low[mask_low]), wave[mask], max(wave_i_high[mask_high])]
            array_i = np.r_[array_i_low[mask_low][0], array[mask], array_i_high[mask_high][-1]]

            return wave_i, array_i


        if interp:
            wave_i, flux_i = interpolate_band_lims(flux, wave, wmin, wmax, step)
            _, flux_err_i = interpolate_band_lims(flux_err, wave, wmin, wmax, step)
        else:
            mask = (wave >= wmin) & (wave <= wmax)
            flux_i = flux[mask]
            flux_err_i = flux_err[mask]
            wave_i = wave[mask]

        r_neg_ln, _, _, _ = check_negflux(flux_i, verb=False)
        
        if bandtype == 'tri':
            bp_i = 1 - np.abs(wave_i - ctr)/win
            bp_i = np.where(bp_i > 0, bp_i, bp_i*0.0)
        elif bandtype == 'sq':
            bp_mask = (wave_i >= ctr - win/2) & (wave_i <= ctr + win/2)
            bp_i = np.where(bp_mask, 1, 0.0)

        if show_plot:
            plt.plot(wave, flux, 'k-', lw=0.7, alpha=0.5)
            plt.plot(wave_i, flux_i, 'k-', lw=1)
            plt.plot(wave_i, flux_i.max()*bp_i, 'b--', lw=0.7, label='bandpass')
            plt.axvline(wmin, c='k', ls=':', lw=0.7, label='win')
            plt.axvline(wmax, c='k', ls=':', lw=0.7)
            plt.xlabel(r"$\lambda$ [$\AA$]")
            plt.ylabel("Norm. Flux")
            plt.legend(loc=1)
            plt.ylim(-1.5*abs(min(flux_i)), 2*max(flux_i))
            plt.xlim(wmin-win/2, wmax+win/2)
            plt.show()

        flux_band = sum(flux_i * bp_i)/win
        flux_band_var= sum((flux_err_i**2 + noise**2) * bp_i**2)/win**2

        return flux_band, np.sqrt(flux_band_var), r_neg_ln
