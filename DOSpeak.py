
import math
import numpy as np
from numpy.ma.core import indices
from scipy import optimize

import plot

verbose = False

def lorentzian(E, y0, A, Gamma, Er)
    return y0 + (A / np.pi) * (Gamma / 2) / ((E - Er) ** 2 + (Gamma / 2) ** 2)


class DOSpeak:
    discontinuity_treatment = "fit"

    def __init__(self, energy, rho, gamma, root: int, peak_E, peak_rho):

        self.root = root
        self.approx_peak_E = peak_E
        self.approx_peak_rho = peak_rho
        self.approx_y0 = min(rho) / 2
        self.energy_array, self.dos_array, self.gamma_array = self.trim(energy, rho, gamma)  # Trimming causes fits to fail, especially if only half the peak is present.
        # self.energy_array, self.dos_array, self.gamma_array = energy, rho, gamma
        self.energy_below, self.energy_above = self.max_energy_window(energy, rho)  # fit_E needs to fall between these
        self.pointwise_energy = self.energy_array[np.argmax(self.dos_array)]
        self.fitted_dos_array = None
        self.fit_E = None
        self.fit_Gamma = None
        self.fit_A = None
        self.fit_y0 = None
        self.fit_gamma = None
        self.ssr = None
        self.rel_ssr_per_point = None
        self.worse = []
        self.trim_left = 0
        self.trim_right = 0
        self.using_pointwise_energy = False
        self.nr_fit_attempts = 0
        self.warning = None


    def trim(self, energy, rho, gamma):
        highest_i = np.argmax(rho)
        left_energy = energy[:highest_i+1]
        left_energy_reverse = left_energy[::-1]
        left_rho = rho[:highest_i+1]
        left_rho_reverse = left_rho[::-1]
        right_energy = energy[highest_i:]
        right_rho = rho[highest_i:]
        trim_left = len(left_rho) - self.trim_half(left_energy_reverse, left_rho_reverse)
        trim_right = len(left_rho) + self.trim_half(right_energy, right_rho)
        self.trim_left = trim_left
        self.trim_right = len(rho) - trim_right
        if verbose:
            if trim_left>0:
                print(f"Root {self.root}, peak {self.approx_peak_E}: {trim_left} points trimmed off left.")
            if trim_right<len(rho):
                print(f"Root {self.root}, peak {self.approx_peak_E}: {len(rho) - trim_right} points trimmed off right.")
        return energy[trim_left:trim_right], rho[trim_left:trim_right], gamma[trim_left:trim_right]


    def trim_half(self, energy, dos):
        deriv_num = (dos[1:] - dos[:-1])
        deriv = deriv_num/(energy[1:] - energy[:-1])
        steepening = True
        prev_delta = 0
        for i, delta in enumerate(deriv):
            if abs(deriv_num[i]) > 0.6*max(dos):
                # print(f"Big jump detected at i={i}! {delta} vs {max(dos)}")
                return i+1
            if steepening and abs(delta) < abs(prev_delta):
                steepening = False  # we're past the steepest point of the peak flank
            if not steepening and abs(delta) > abs(prev_delta):
                # print(f"i={i}, re-steepening cut-off")
                return i  # getting steeper again; cut off here.
            prev_delta = delta
        return len(dos)  # no trimming condition found; use whole arrays


    def get_smooth_lorentzian_curve(self, x_array):
        return lorentzian(x_array, self.fit_y0, self.fit_A, self.fit_Gamma, self.fit_E)

    def print_fitted_parameters(self):
        print(f"Root #{self.root}, E = {self.fit_E}, Gamma = {self.fit_Gamma}, A = {self.fit_A}, y0 = {self.fit_y0}")

    def estimate_gamma(self):
        half_max = 0.5 * self.approx_peak_rho
        max_index = np.argmax(self.dos_array)
        dos_left = self.dos_array[:max_index]
        dos_right = self.dos_array[max_index:]
        # indices = np.where(self.dos_array > half_max)[0]
        indices_left = np.where(dos_left > half_max)[0]
        indices_right = np.where(dos_right > half_max)[0]
        if len(indices_left)+len(indices_right) < 2:
            return self.energy_array[1] - self.energy_array[0]  # not enough data points, rough guess
        # return abs(self.energy_array[indices[-1]] - self.energy_array[indices[0]])
        left_width = abs(self.approx_peak_E - self.energy_array[indices_left[0]]) if len(indices_left) else 0
        right_width = abs(self.approx_peak_E - (self.energy_array[indices_right[-1] + max_index])) if len(indices_right) else 0
        return 2*max(left_width, right_width)  # Deal with half-peak cases
        # return abs(self.energy_array[indices[-1]] - self.energy_array[indices[0]])

    def initial_guesses(self):
        """
        Provide initial guesses for the fitting parameters.

        Returns:
        list: Initial guesses [y0, A, Gamma, Er] for the Lorentzian fit.
        """
        y0 = self.approx_y0
        Gamma = self.estimate_gamma()
        A = self.approx_peak_rho * np.pi * Gamma / 2
        Er = self.approx_peak_E
        mid_index = int(len(self.dos_array)/2)
        y0 += self.dos_array[mid_index] - lorentzian(self.energy_array[mid_index], y0, A, Gamma, Er)

        return [y0, A, Gamma, Er]


    def fit_lorentzian(self):
        """
        Fit the Lorentzian model to the DOS data.

        Returns:
        list: Optimized parameters [y0, A, Gamma, Er] or None if the fit failed.
        """
        if self.root == 14:
            y0, A, Gamma, Er = self.initial_guesses()
            guessed_peak = DOSpeak(self.energy_array, self.dos_array, self.gamma_array, self.root, self.approx_peak_E, self.approx_peak_rho)
            guessed_peak.fit_E = Er
            guessed_peak.fit_A = A
            guessed_peak.fit_y0 = y0
            guessed_peak.fit_Gamma = Gamma
            plot.peak_fit(guessed_peak, f"/home/giogina/Desktop/DOSmax/stb_plot_pf/{Er}.png")
        self.nr_fit_attempts += 1
        try:
            guesses = self.initial_guesses()
            if self.nr_fit_attempts == 1:
                cv = optimize.curve_fit(lorentzian, self.energy_array, self.dos_array, p0=guesses)
            else:
                cv = optimize.curve_fit(lorentzian, np.insert(self.energy_array, 0, self.approx_peak_E-5*guesses[2]),
                                        np.insert(self.dos_array, 0, guesses[0]), p0=guesses)
        except Exception as e:
            if verbose:
                print(f"Root {self.root}: Fit failed for peak at E={self.approx_peak_E}, rho={self.approx_peak_rho}, on {len(self.energy_array)} available data points.")
            if 2 < self.nr_fit_attempts < 9 and len(self.dos_array) > 20:
                if verbose: print("Retrying...")  # todo: also retry above?
                self.approx_peak_E -= 0.0001  # todo: less arbitrary steps
                self.approx_peak_rho -= 0.0001 * (self.dos_array[2] - self.dos_array[1])/(self.energy_array[2] - self.energy_array[1])
            if self.nr_fit_attempts < 9 and len(self.dos_array) > 20:
                return self.fit_lorentzian()
            else:
                self.using_pointwise_energy = True
                return None
        popt = cv[0]
        pcov = cv[1]
        if np.any(np.isinf(pcov)):
            if verbose:
                print(f"Root {self.root}: Fit did not converge for peak at E={self.approx_peak_E}, rho={self.approx_peak_rho}, on {len(self.energy_array)} available data points.")
            return None
        else:
            self.fitted_dos_array = np.array([lorentzian(e, *popt) for e in self.energy_array])
            self.ssr = np.sum((self.dos_array - self.fitted_dos_array) ** 2)
            self.rel_ssr_per_point = math.sqrt(self.ssr) / (len(self.dos_array)) / self.approx_peak_rho
            self.fit_E = popt[3]
            self.fit_Gamma = popt[2]
            self.fit_A = popt[1]
            self.fit_y0 = popt[0]
            self.fit_gamma = np.interp(self.fit_E, self.energy_array, self.gamma_array)

            if verbose:
                print(f"Root {self.root}: Er = {int(popt[3] * 10000) / 10000}, Gamma = {int(popt[2] * 1000000) / 1000000}, A = {int(popt[1] * 10000) / 10000}, y0 = {int(popt[0] * 1000) / 1000};     Relative SSR per data point: {self.rel_ssr_per_point}")

            self.check_fit()

            return popt

    def max_energy_window(self, energy, rho):
        """
        Determines the energy window surrounding the maximum DOS value.

        The function finds the maximum value in the rho array and returns the energies
        of the points directly to the left and right of this maximum. These surrounding
        energy values define the range (window) within which the fitted energy (fit_E)
        should fall to be considered "close enough" to the maximum DOS.

        Parameters:
        energy (numpy.ndarray): The array of energy values.
        rho (numpy.ndarray): The array of DOS values corresponding to the energies.

        Returns:
        tuple: A pair (e_below, e_above) representing the energies directly to the left
               and right of the maximum DOS value. If the maximum is at the first or last
               position, e_below or e_above will be None, respectively.
        """
        # Find the index of the maximum DOS value
        max_index = np.argmax(rho)

        # Handle edge cases where the maximum is at the first or last position
        if max_index == 0:
            e_below = None
            e_above = energy[max_index + 1]
        elif max_index == len(energy) - 1:
            e_below = energy[max_index - 1]
            e_above = None
        else:
            e_below = energy[max_index - 1]
            e_above = energy[max_index + 1]
        if e_below is None or e_above is None:
            print(f"Warning: Root {self.root}, peak {self.approx_peak_E}: Maximum is at edge.")

        return e_below, e_above

    def check_fit(self):
        max_index = np.argmax(self.dos_array)
        if max_index == 0:
            self.warning = "Only right side of peak available to fit!"
            if self.discontinuity_treatment != "fit":
                self.warning += " Using pointwise maximum."
                self.using_pointwise_energy = True
            return True
        if max_index == len(self.dos_array)-1:
            self.warning = "Only left side of peak available to fit!"
            if self.discontinuity_treatment != "fit":
                self.warning += " Using pointwise maximum."
                self.using_pointwise_energy = True
            return True
        if self.energy_below is not None and self.energy_below > self.fit_E:
            self.warning = "Fit maximum not in range of pointwise maximum. Likely bad fit."
            return False
        if self.energy_above is not None and self.energy_above < self.fit_E:
            self.warning = "Fit maximum not in range of pointwise maximum. Likely bad fit."
            return False
        return True

    def energy(self):
        if DOSpeak.discontinuity_treatment == "fit":
            return self.fit_E
        else:
            return self.pointwise_energy if self.using_pointwise_energy else self.fit_E