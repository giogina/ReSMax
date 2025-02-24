import math
import numpy as np
import warnings
from scipy import optimize
from scipy.signal import argrelextrema

verbose = False


def lorentzian(E, y0, A, Gamma, Er):
    """
    Calculate the Lorentzian function.

    Parameters:
    E (float): The energy value.
    y0 (float): The baseline offset.
    A (float): The amplitude of the peak.
    Gamma (float): The full width at half maximum (FWHM).
    Er (float): The resonance energy (peak center).

    Returns:
    float: The Lorentzian function value at energy E.
    """
    return y0 + (A / np.pi) * (Gamma / 2) / ((E - Er) ** 2 + (Gamma / 2) ** 2)


class DOSpeak:
    """
    Class representing a Density of States (DOS) peak and handling the fitting of the peak
    to a Lorentzian model.
    """
    discontinuity_treatment = "fit"

    def __init__(self, energy, rho, gamma, root: int):
        """
        Initialize the DOSpeak object.

        Parameters:
        energy (array-like): Array of energy values.
        rho (array-like): Array of density of states values.
        root (int): The root identifier for this DOS peak.
        peak_E (float): The approximate energy of the peak.
        peak_rho (float): The approximate density of states at the peak.
        """

        self.root = root
        self.is_descending = (min(rho) < 0)
        peak_i = rho.argmax() if not self.is_descending else argrelextrema(energy, np.less)[0][0] if len(argrelextrema(energy, np.less)[0]) > 0 else 0
        self.approx_peak_E = float(energy[peak_i])
        self.approx_peak_rho = float(rho[peak_i])
        self.approx_y0 = min(rho) / 2
        self.approx_Gamma = None
        self.approx_A = None
        self.energy_array, self.dos_array, self.gamma_array = self.trim(energy, rho, gamma)
        self.pointwise_energy = self.energy_array[np.argmax(self.dos_array)]
        # self.slope_energies = []  # temp
        self.deriv_ratio = self.compute_dos_derivative_ratio()
        self.fitted_dos_array = None
        self.fit_E = None
        self.fit_Gamma = None
        self.fit_A = None
        self.fit_y0 = None
        self.fit_gamma = float(gamma[peak_i])
        self.ssr = None
        self.rel_ssr_per_point = None
        self.trim_left = 0
        self.trim_right = 0
        self.using_pointwise_energy = False
        self.nr_fit_attempts = 0
        self.energy_below, self.energy_above = None, None  # fit_E needs to fall between these
        self.warning = None

    def trim(self, energy, rho, gamma):
        if self.is_descending:
            return energy, rho, gamma  # this is a descending piece
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
        dos_max = np.max(dos)
        deriv_num = (dos[1:] - dos[:-1])
        big_jump_indices = np.where(np.abs(deriv_num) > 0.6 * dos_max)[0]
        if big_jump_indices.size > 0:
            return big_jump_indices[0] + 1  # Return first occurrence of a jump
        return len(dos)  # no trimming condition found; use whole arrays

    def get_smooth_lorentzian_curve(self, x_array):
        """
        Generate a smooth Lorentzian curve based on fitted parameters.

        Parameters:
        x_array (array-like): Array of x values over which to calculate the Lorentzian curve.

        Returns:
        array-like: Smooth Lorentzian curve values corresponding to x_array.
        """
        return lorentzian(x_array, self.fit_y0, self.fit_A, self.fit_Gamma, self.fit_E)

    def get_smooth_guess_lorentzian_curve(self, x_array):  # temp: debug plotting
        return lorentzian(x_array, self.approx_y0, self.approx_A, self.approx_Gamma, self.approx_peak_E)

    def print_fitted_parameters(self):
        """
        Print the fitted parameters for the DOS peak.
        """
        print(f"Root #{self.root}, E = {self.fit_E}, Gamma = {self.fit_Gamma}, A = {self.fit_A}, y0 = {self.fit_y0}")

    def distance_in_Gamma(self, peak):
        """
        Get the distance to another peak's maximum as multiple of own Gamma
        """
        if self.fit_E is None or peak.fit_E is None:
            return None
        return np.abs(self.fit_E - peak.fit_E) / self.fit_Gamma

    def estimate_gamma(self):
        """
        Estimate the Gamma parameter (width of the peak) based on the DOS data.

        Returns:
        float: Estimated Gamma value.
        """
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

    def make_initial_guesses(self):
        """
        Provide initial guesses for the fitting parameters; save in self.approx_* class instance variables.
        """

        self.approx_Gamma = self.estimate_gamma()
        self.approx_A = self.approx_peak_rho * np.pi * self.approx_Gamma / 2
        mid_index = int(len(self.dos_array)/2)
        self.approx_y0 += self.dos_array[mid_index] - lorentzian(self.energy_array[mid_index], self.approx_y0, self.approx_A, self.approx_Gamma, self.approx_peak_E)


    def fit_lorentzian(self):
        """
        Fit the Lorentzian model to the DOS data.

        Returns:
        list: Optimized parameters [y0, A, Gamma, Er] or None if the fit failed.
        """
        if self.is_descending:
            self.rel_ssr_per_point = 10**6 + self.approx_peak_E  # large value for fit quality
            return None
        if self.deriv_ratio < 0.3:
            if verbose:
                print(f"Root {self.root}: Peak at E={self.approx_peak_E}, rho={self.approx_peak_rho} is too asymmetrical to bother fitting.")
            return None
        self.nr_fit_attempts += 1
        try:
            self.make_initial_guesses()
            self.trim_arrays_to_n_Gammas_around_max()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", optimize.OptimizeWarning) #(E, y0, A, Gamma, Er)
                if self.nr_fit_attempts == 1:
                    cv = optimize.curve_fit(lorentzian, self.energy_array, self.dos_array, p0=[self.approx_y0, self.approx_A, self.approx_Gamma, self.approx_peak_E])
                else:
                    energy_prepended = np.concatenate(([self.approx_peak_E - 5 * self.approx_Gamma], self.energy_array))
                    dos_prepended = np.concatenate(([self.approx_y0], self.dos_array))
                    cv = optimize.curve_fit(lorentzian, energy_prepended, dos_prepended, p0=[self.approx_y0, self.approx_A, self.approx_Gamma, self.approx_peak_E])

        except Exception as e:
            if verbose:
                print(f"Root {self.root}: Fit failed for peak at E={self.approx_peak_E}, rho={self.approx_peak_rho}, on {len(self.energy_array)} available data points.")
            if 2 < self.nr_fit_attempts < 3 and len(self.dos_array) > 20:
                if verbose: print("Retrying...")
                self.approx_peak_E -= 0.0001
                self.approx_peak_rho -= 0.0001 * (self.dos_array[2] - self.dos_array[1])/(self.energy_array[2] - self.energy_array[1])
            if self.nr_fit_attempts < 9 and len(self.dos_array) > 20:
                return self.fit_lorentzian()
            else:
                self.using_pointwise_energy = True
                return None
        popt = cv[0]
        pcov = cv[1]
        # if self.approx_peak_E < -0.55:  # looks like the guesses are pretty good
        #     print((popt-guesses)/popt, [float(g) for g in guesses], [float(g) for g in popt])  # [y0, A, Gamma, Er]
        if np.any(np.isinf(pcov)):
            if verbose:
                print(f"Root {self.root}: Fit did not converge for peak at E={self.approx_peak_E}, rho={self.approx_peak_rho}, on {len(self.energy_array)} available data points.")
            return None
        else:
            self.fitted_dos_array = lorentzian(self.energy_array, *popt)
            self.ssr = np.sum((self.dos_array - self.fitted_dos_array) ** 2)
            self.rel_ssr_per_point = math.sqrt(self.ssr) / (len(self.dos_array)) / self.approx_peak_rho**2  # / (0.01+self.deriv_ratio)
            self.fit_E = popt[3]
            self.fit_Gamma = popt[2]
            self.fit_A = popt[1]
            self.fit_y0 = popt[0]
            self.fit_gamma = np.interp(self.fit_E, self.energy_array, self.gamma_array)

            if verbose:
                print(f"Root {self.root}: Er = {int(popt[3] * 10000) / 10000}, Gamma = {int(popt[2] * 1000000) / 1000000}, A = {int(popt[1] * 10000) / 10000}, y0 = {int(popt[0] * 1000) / 1000};     Relative SSR per data point: {self.rel_ssr_per_point}")

            self.check_fit()
            # print(f"E: {abs((self.fit_E - self.approx_peak_E)/self.fit_E)*100:.3f}%, Gamma: {abs((self.fit_Gamma - guesses[2])/self.fit_Gamma)*100:.3f}%, {self.rel_ssr_per_point:.6f}, {self.warning}")
            return popt


    def trim_arrays_to_n_Gammas_around_max(self, multiple = 10):
        """
        Trim the energy, DOS, and gamma arrays to fit within [fit_E - multiple*fit_Gamma, fit_E + multiple*fit_Gamma].
        """

        energy = self.approx_peak_E
        gamma = self.approx_Gamma
        if energy is None or gamma is None:
            return

        E_min = energy - multiple * gamma
        E_max = energy + multiple * gamma

        # Find indices within the trim range
        energy_mask = (self.energy_array >= E_min) & (self.energy_array <= E_max)
        valid_indices = np.where(energy_mask)[0]

        if len(valid_indices) == 0:
            return

        trim_front = valid_indices[0]  # How many points are removed from the front
        trim_back = len(self.energy_array) - valid_indices[-1] - 1  # From the back

        self.energy_array = self.energy_array[valid_indices]
        self.gamma_array = self.gamma_array[valid_indices]

        # Trim DOS array consistently (it has one fewer point)
        if len(self.dos_array) > trim_front + trim_back:
            self.dos_array = self.dos_array[trim_front:-trim_back] if trim_back > 0 else self.dos_array[trim_front:]
        else:
            return

        if verbose:
            print(f"Root {self.root}: Trimmed arrays to range {E_min:.6f} - {E_max:.6f}, removed {trim_front} from front, {trim_back} from back.")


    def compute_dos_derivative_ratio(self):
        """
        Compute the ratio of the smallest to largest absolute value of the
        minimum and maximum derivatives of the DOS vs. energy.

        Returns:
        float: The ratio of |min_derivative| to |max_derivative| for DOS vs. energy.
        """
        if self.is_descending:
            return 0
        max_index = np.argmax(self.dos_array)
        if max_index < 3 or max_index > len(self.dos_array) - 4:
            return 0  # If too few points are available on one side, return 0.

        derivatives = np.gradient(self.dos_array, self.energy_array)
        min_derivative = abs(np.min(derivatives))
        max_derivative = abs(np.max(derivatives))

        # max_slope_index = np.argmax(derivatives)
        # max_slope = derivatives[max_slope_index]
        # if len(derivatives[:max_slope_index]) and len(derivatives[max_slope_index:]):
        #     half_slope_1 = np.argmin(np.abs(derivatives[max_slope_index:] - max_slope / 2))+max_slope_index
        #     half_slope_2 = np.argmin(np.abs(derivatives[:max_slope_index] - max_slope / 2))
        #     # Em = self.energy_array[max_slope_index]
        #     Ehs1 = self.energy_array[half_slope_1]
        #     Ehs2 = self.energy_array[half_slope_2]
        #     Er_by_max = self.energy_array[max_index]
        #     Er_by_outer_half_slope = 1.424839036*Ehs1-.4248390359*Em
        #     Er_by_inner_half_slope = -.7136756879*Ehs2+1.713675688*Em
        #     W_by_outer = 4.935787206*(Ehs1-Em)
        #     W_by_inner = 2.472245103*(Em-Ehs2)
        #     print(f"Er guesses: {Er_by_max}, {Er_by_outer_half_slope}, {Er_by_inner_half_slope};  Er guess diffs: {(Er_by_outer_half_slope-Er_by_max)/(W_by_outer/2+W_by_inner/2)*100:.1f}%, {(Er_by_inner_half_slope-Er_by_max)/(W_by_outer/2+W_by_inner/2)*100:.1f}%, {(Er_by_outer_half_slope-Er_by_inner_half_slope)/(W_by_outer/2+W_by_inner/2)*100:.1f}%")
        #     print(f"Gamma guesses: {W_by_outer}, {W_by_inner}")
        #     self.slope_energies = [max_slope_index, half_slope_1, half_slope_2, self.energy_array[max_index], 1.424839036*Ehs1-.4248390359*Em,-.7136756879*Ehs2+1.713675688*Em, f"Er guess diffs: {(Er_by_outer_half_slope-Er_by_max)/(W_by_outer/2+W_by_inner/2)*100:.1f}%, {(Er_by_inner_half_slope-Er_by_max)/(W_by_outer/2+W_by_inner/2)*100:.1f}%, {(Er_by_outer_half_slope-Er_by_inner_half_slope)/(W_by_outer/2+W_by_inner/2)*100:.1f}%"]

        if min(min_derivative, max_derivative) != 0:
            ratio = min(min_derivative, max_derivative) / max(min_derivative, max_derivative)
        else:
            ratio = 0

        # print(f"Min derivative = {min_derivative:.6f}")
        # print(f"Max derivative = {max_derivative:.6f}")
        # print(f"Ratio (|min|/|max|) = {ratio:.6f}")

        return ratio

    def max_energy_window(self):
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

        if self.fit_Gamma is not None:
            e_below = self.pointwise_energy - self.fit_Gamma*0.2
            e_above = self.pointwise_energy + self.fit_Gamma*0.2
        else:
            e_below = None
            e_above = None

        # max_index = np.argmax(self.dos_array)
        # if max_index == 0:
        #     e_below = None
        #     e_above = self.energy_array[max_index + 1]
        # elif max_index == len(self.energy_array) - 1:
        #     e_below = self.energy_array[max_index - 1]
        #     e_above = None
        # else:
        #     e_below = self.energy_array[max_index - 1]
        #     e_above = self.energy_array[max_index + 1]

        # if self.fit_Gamma is not None:
        #     print(e_below, e_above, self.pointwise_energy - self.fit_Gamma*0.01, self.pointwise_energy + self.fit_Gamma*0.01)
        # if e_below is None or e_above is None:
        #     print(f"Warning: Root {self.root}, peak {self.approx_peak_E}: Maximum is at edge.")

        return e_below, e_above

    def check_fit(self):
        self.energy_below, self.energy_above = self.max_energy_window()
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
        if self.is_descending:
            return self.approx_peak_E
        return self.fit_E
            
