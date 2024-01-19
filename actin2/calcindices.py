import sys, os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from scipy.interpolate import interp1d


from .indtable import IndTable
from .spectrographs._spec_tools import printif


class CalcIndices:
    """Calculates activity indices.

    Args:
        spectrum (dict): Dictionary containing wavelength ``wave`` and flux ``flux`` data.
        headers (dict): Dictionary containing useful fits headers information.
        indices (list, str): List of indices ``ìnd_id`` to be extracted.
        step (float): Step used in interpolation. Must have ``interp=True``
        table_df (pd.DataFrame, None): Indices table to be used as input. If ``None`` the default ``actin_table.csv`` will be used.
        plot_lines (bool): Diagnostic plot.
        full_output (bool): if ``True`` the output includes the fluxes and errors in of each line.
        interp (bool): Use interpolation to reeduce the effect of finite resolution of thhe spectrographs. The interpolation step is ``step`` value.
        verb (bool): Turn verbose on/off.

    Attributes:
        indices (dict): Dictionary with indices data.
    """

    def __init__(self, spectrum, headers, indices, step=1e-5, table_df=None, plot_lines=False, full_output=False, interp=True, verb=False):

        if "flux" not in spectrum:
            raise KeyError(f"{__class__.__name__} ERROR: There is no 'flux' keyword in the 'spectrum' dictionary")

        if "wave" not in spectrum:
            raise KeyError(f"{__class__.__name__} ERROR: There is no 'wave' keyword in the 'spectrum' dictionary")

        # If indices is given as a string enforce type as list:
        if not isinstance(indices, (list, np.ndarray, tuple)):
            indices = [indices]

        if table_df is None:
            ind_table = IndTable().table
        else:
            if not isinstance(table_df, pd.DataFrame):
                raise RuntimeError(f"{__class__.__name__} ERROR: 'table_df' is not a pandas DataFrame")
            ind_table = table_df

        # Confirm that index and required lines (at least one activity and reference line per index) are present in the table:
        for index in indices:
            if index not in np.unique(ind_table.ind_id.values):
                raise RuntimeError(f"{__class__.__name__} ERROR: Index {index} not available in 'ind_table'. Available indices are: {np.unique(ind_table.ind_id.values)}")
            if "R1" not in ind_table[ind_table.ind_id==index].ind_var.values:
                raise RuntimeError(f"{__class__.__name__} ERROR: Index {index} does not have reference lines ('R1' under 'ind_var' column) in 'ind_table'. Available line variables for index are: {np.unique(ind_table[ind_table.ind_id==index].ind_var.values)}")
            if "L1" not in ind_table[ind_table.ind_id==index].ind_var.values:
                raise RuntimeError(f"{__class__.__name__} ERROR: Index {index} does not have activity lines ('L1' under 'ind_var' column) in 'ind_table'. Available line variables for index are: {np.unique(ind_table[ind_table.ind_id==index].ind_var.values)}")
                

        flux = spectrum['flux']
        wave = spectrum['wave']


        if "flux_err" in spectrum:
            flux_err = spectrum['flux_err'] # flux noise from spectrograph module, depends on spectrograph
        else:
            flux_err = np.sqrt(abs(flux)) # photon noise


        # Usefull if a 'noise' key with additional noise comes from the spectrograph 'header' dictionary with additional noise value to be added in quadrature to the index errors:
        if 'noise' in headers:
            noise = headers['noise']
        else:
            noise = 0.0


        data = dict()

        for index in indices:
            index_dict = self.calc_index(wave, flux, flux_err, noise, index, step, ind_table, plot_lines, full_output, interp, verb)

            data.update(index_dict)

        self.indices = data


    def calc_index(self, wave, flux, flux_err, noise, index, step, ind_table, plot_lines, full_output, interp, verb):

        # get index from ind_table:
        df = ind_table[ind_table.ind_id == index]
        ln_ids = df.ln_id.values


        F_num = np.zeros(len(ln_ids))
        F_num_err = np.zeros(len(ln_ids))
        F_denum = np.zeros(len(ln_ids))
        F_denum_err = np.zeros(len(ln_ids))
        Rneg_lines = np.zeros(len(ln_ids))
        win_list = np.zeros(len(ln_ids))

        lines = dict()

        for i, line in enumerate(ln_ids):
            mask = (df.ln_id == line)

            ctr = df.ln_ctr[mask].values[0]
            win = df.ln_win[mask].values[0]
            bandtype = df.bandtype[mask].values[0]
            const = df.ln_c[mask].values[0]
            ind_var = df.ind_var[mask].values[0]


            if len(wave.shape) == 2:
                order = self._spec_order(wave, flux, flux_err, ctr, win, bandtype, verb=verb)
                w = wave[order]
                f = flux[order]
                f_err = flux_err[order]
            if len(wave.shape) == 1:
                w = self._check_wave_range(wave, index, ind_table, verb=verb)
                f = flux; f_err = flux_err


            if w is None or f is None: # the exception was taken care off above
                return dict()

            flux_band, flux_band_err, Rneg_ln, npix = self.calc_flux_band(w, f, f_err, noise, ctr, win, bandtype, step, show_plot=plot_lines, interp=interp, verb=verb)

            if full_output:
                lines[line + '_F'] = flux_band
                lines[line + '_F_err'] = flux_band_err
                lines[line + '_Rneg'] = Rneg_ln
                lines[line + '_npix'] = npix
                if len(wave.shape) == 2:
                    lines[line + '_order'] = order
                

            Rneg_lines[i] = Rneg_ln

            if ind_var.startswith("L"):
                F_num[i] = const * flux_band
                F_num_err[i] = flux_band_err

            elif ind_var.startswith("R"):
                F_denum[i] = const * flux_band
                F_denum_err[i] = flux_band_err

            win_list[i] = win


        ind = np.sum(F_num)/np.sum(F_denum)

        ind_err = np.sqrt(np.sum(F_num_err**2) + ind**2 * np.sum(F_denum_err**2)) / (np.sum(F_denum))

        # Ratio of negative flux to total flux:
        Rneg = np.mean(Rneg_lines)

        data = {}
        data[index] = ind
        data[index + '_err'] = ind_err
        data[index + '_Rneg'] = Rneg

        if full_output:
            data.update(lines)

        return data


    def calc_flux_band(self, wave, flux, flux_err, noise, ctr, win, bandtype, step, show_plot=False, interp=True, verb=False):

        if bandtype == "tri":
            wmin = ctr - win; wmax = ctr + win
        if bandtype == "sq":
            wmin = ctr - win/2; wmax = ctr + win/2

        def _interp_band_lims(array, wave, wmin, wmax, step):
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
            # wavelength right before and after the window limits:
            wave_int_low = (wave[wave < wmin][-1], wave[wave >= wmin][0])
            wave_int_high = (wave[wave <= wmax][-1], wave[wave > wmax][0])

            reso_low = np.array(wave_int_low).ptp()
            reso_high = np.array(wave_int_high).ptp()

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

            mask_low = (wave_i_low >= wmin)
            mask_high = (wave_i_high <= wmax)

            frac_low = wave_i_low[mask_low].ptp()/reso_low
            frac_high = wave_i_high[mask_high].ptp()/reso_high

            mask = (wave >= wmin) & (wave <= wmax)

            wave_i = np.r_[min(wave_i_low[mask_low]), wave[mask], max(wave_i_high[mask_high])]
            array_i = np.r_[array_i_low[mask_low][0] * frac_low, array[mask], array_i_high[mask_high][-1] * frac_high]

            return wave_i, array_i, frac_low, frac_high


        if interp:
            wave_i, flux_i, frac_low, frac_high = _interp_band_lims(flux, wave, wmin, wmax, step)
            _, flux_err_i, _, _ = _interp_band_lims(flux_err, wave, wmin, wmax, step)
            npix = len(flux_i) -2 + frac_low + frac_high
        else:
            mask = (wave >= wmin) & (wave <= wmax)
            flux_i = flux[mask]
            flux_err_i = flux_err[mask]
            wave_i = wave[mask]
            npix = len(flux_i)

        Rneg_ln = self._calc_Rneg_ln(flux_i)

        
        if bandtype == 'tri':
            bp_i = 1 - np.abs(wave_i - ctr)/win
            bp_i = np.where(bp_i > 0, bp_i, bp_i*0.0)
        elif bandtype == 'sq':
            bp_mask = (wave_i >= ctr - win/2) & (wave_i <= ctr + win/2)
            bp_i = np.where(bp_mask, 1, 0.0)

        flux_band = sum(flux_i * bp_i)/win
        flux_band_var = sum((flux_err_i**2 + noise**2) * bp_i**2)/win**2

        # plt.title(f"npix = {npix}")
        # plt.plot(wave, flux, 'k.-')
        # plt.plot(wave_i, flux_i, 'b.-')
        # plt.plot(wave_i, bp_i, 'g.-')
        # plt.axvline(wmin, color='r')
        # plt.axvline(wmax, color='r')
        # plt.show()

        if show_plot:
            if bandtype == 'tri':
                bp = 1 - np.abs(wave - ctr)/win
                bp = np.where(bp > 0, bp, bp*0.0)
            elif bandtype == 'sq':
                bp_mask = (wave >= ctr - win/2) & (wave <= ctr + win/2)
                bp = np.where(bp_mask, 1, 0.0)

            plt.plot(wave, flux, 'k-', lw=0.7, alpha=0.5)
            plt.plot(wave_i, flux_i, 'k-', lw=1)
            plt.plot(wave, 1.2*flux_i.max()*bp, 'b--', lw=0.7, label='bandpass')

            plt.xlabel(r"$\lambda$ [$\mathrm{\AA{}}$]")
            plt.ylabel("Norm. Flux")
            plt.legend(loc=1)
            plt.ylim(-0.5*abs(min(flux_i)), 1.5*max(flux_i))
            plt.xlim(wmin-win/2, wmax+win/2)
            plt.show()

        return flux_band, np.sqrt(flux_band_var), Rneg_ln, npix #!


    def _check_wave_range(self, wave, index, table_df, verb=False):
        ln_ids = table_df[table_df.ind_id==index].ln_id.values

        for line in ln_ids:
            ctr = table_df[table_df.ln_id==line].ln_ctr.values[0]
            win = table_df[table_df.ln_id==line].ln_win.values[0]
            bt = table_df[table_df.ln_id==line].bandtype.values[0]

            if bt == 'sq':
                win = win/2

            if wave[0] >= ctr - win or wave[-1] <= ctr + win:
                printif(f"{sys._getframe().f_code.co_name} ERROR: bandpass outside wavelength range", verb=verb)
                return None
            else:
                return wave


    @staticmethod
    def _spec_order(wave_2d, flux_2d, flux_2d_err, ln_ctr, ln_win, bandtype, type='dist', show_orders=False, verb=False):

        if flux_2d_err is None:
            flux_2d_err = np.zeros_like(flux_2d)

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


        if not orders:
            printif(f"{sys._getframe().f_code.co_name} ERROR: bandpass outside wavelength range", verb=verb)
            return None, None, None

        printif(f"available orders for line at {ln_ctr:.3f}: {orders}", verb)

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
        if type == 'first':
            order = orders[0]
            printif(f"selected: {order}", verb)
        if type == 'second':
            order = orders[1]
            printif(f"selected: {order}", verb)

        #return wave_2d[order], flux_2d[order], flux_2d_err[order]
        return order



    def _calc_Rneg_ln(self, flux):
        r"""Negative flux ratio.

        .. Math::

            R_\mathrm{neg} = \frac{\sum(F_\mathrm{neg})} {\sum(F_\mathrm{tot})}


        Args:
            flux (array): Flux inside a bandpass.

        Returns:
            float: Negative flux ratio for a given spectral line [0, 1]
        """
        flux = np.asarray(flux)
        return abs(np.sum(flux[flux < 0]))/np.sum(abs(flux))
