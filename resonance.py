import numpy as np
from scipy import signal

from DOSpeak import DOSpeak

verbose = False

class Resonance:
    """
    Class to represent and manage a resonance and its associated peaks.
    """
    resonances = []

    def __init__(self, peak):
        """
        Initialize the Resonance object.

        Parameters:
        peak (DOSpeak or list): A DOSpeak object or a list of DOSpeak objects associated with this resonance.
        """
        if type(peak) is list:
            self.peaks = peak
            self.best_fit = None
            self.energy = None
        else:
            self.peaks = [peak]
            self.energy = peak.energy()
            self.best_fit = peak
        self.select_best_fit()
        self.index = len(self.resonances)
        self.resonances.append(self)
        self.threshold = None
        self.manual_peak_selection = False  # to show manually selected resonances even if they don't fit the criteria

    def select_best_fit(self):
        fit_qualities = [p.rel_ssr_per_point for p in self.peaks]
        self.best_fit = self.peaks[fit_qualities.index(min(fit_qualities))]
        self.energy = self.best_fit.energy()

    def categorize_by_thresholds(self, thresholds):
        """
        Categorize the resonance by given energy thresholds.

        Parameters:
        thresholds (list): List of threshold energy values.
        """
        above = [t for t in thresholds if t > self.energy]
        if len(above):
            self.threshold = min(above)

    def should_be_shown(self):  # Show this resonance in lists&overview plot?
        if self.best_fit is not None and self.manual_peak_selection:  # always include manually chosen resonances
            return None if self.best_fit.is_descending else True
        if self.best_fit is None or len(self.peaks) < 3:  # No good fit / manually deselected
            return False
        if not self.best_fit.is_descending:  # at least one properly fit section
            return True
        if self.best_fit.is_descending and len(self.peaks) >= 3:  # Several descending peaks; show with warning
            return None
        return False

#
# def energy_clustering(data, e_min, e_max, resolution):
#     """
#     Detect energy clustering by creating a weighted scan across the energy range.
#
#     Parameters:
#     data (dict): Parsed data dictionary containing "gamma", "rho_{root}" arrays, and others.
#     e_min (float): Minimum energy value for the scan.
#     e_max (float): Maximum energy value for the scan.
#     resolution (float): Energy resolution (step size) for the scan.
#
#     Returns:
#     tuple: (energy_grid, clustering_array) where:
#         - energy_grid: The array of energy values used for the scan.
#         - clustering_array: The weighted clustering values for each energy bin.
#     """
#     # Create an energy grid spanning the range [e_min, e_max] with the given resolution
#     energy_grid = np.arange(e_min, e_max, resolution)
#     clustering_array = np.zeros_like(energy_grid)
#
#     kernel_half_width = 9  # Number of bins a point influences
#     kernel = np.exp(-0.5 * (np.arange(-kernel_half_width, kernel_half_width+1) / kernel_half_width * 3) ** 2)  # Gaussian kernel
#     kernel /= kernel.sum()
#
#     for key in data.keys():
#         if isinstance(key, str) and key.startswith("rho_"):  # Identify DOS arrays
#             rho = data[key]
#             energy = data[int(key[4:])]  # Corresponding energy array for the root
#             log_dos = np.log10(rho+1)
#
#             # Loop over all data points and distribute their weights to the energy grid
#             for e, log_rho in zip(energy[1:-1], log_dos):  # Skip endpoints
#                 # Find the closest grid index
#                 idx = int((e - e_min) / resolution)
#                 if idx < 0 or idx >= len(energy_grid):
#                     continue  # Skip points outside the range
#
#                 # Distribute weights across kernel_width bins
#                 for offset, weight in enumerate(kernel):
#                     bin_idx = idx + offset - kernel_half_width
#                     if 0 <= bin_idx < len(clustering_array):
#                         clustering_array[bin_idx] += log_rho * weight
#
#     peak_indices, _ = signal.find_peaks(
#         clustering_array,
#         # width=2 * resolution / (energy_grid[1] - energy_grid[0]),  # Min peak width in bins
#         height=1
#     )
#     # valid_peak_indices = filter_peaks_by_relative_prominence(clustering_array, peak_indices, prominence_threshold=0.3)
#     valid_peak_indices = peak_indices
#
#     peak_energies = energy_grid[valid_peak_indices] if len(valid_peak_indices) else []
#     return energy_grid, clustering_array, peak_energies

def find_resonances(peaks_by_root):
    """
    Find resonances by analyzing the peaks from all roots.

    Parameters:
    peaks_by_root (dict): Dictionary of peaks organized by root.
    """
    all_peaks = [peak for peaks in peaks_by_root.values() for peak in peaks]
    all_peaks.sort(key=lambda p: p.energy())

    waiting_peaks = []
    for i, peak in enumerate(all_peaks):
        waiting_peaks.append(peak)
        if peak.root in [p.root for p in waiting_peaks[:-1]]:  # repetition found; need to draw a line somewhere between the waiting peaks
            if verbose:
                print(f"Repetition: {peak.root} in {[p.root for p in waiting_peaks[:-1]]}")
            rep_index = [p.root for p in waiting_peaks[:-1]].index(peak.root)
            dists = [waiting_peaks[i + 1].energy() - waiting_peaks[i].energy() if i >= rep_index else 0 for i in range(0, len(waiting_peaks) - 1)]  # 0 before rep index ensures cut happens after.
            max_jump = dists.index(max(dists))
            res = Resonance(waiting_peaks[:max_jump + 1])
            if verbose:
                print("~~~~~", res.energy)
                for p in waiting_peaks[:max_jump + 1]:
                    print(f"    {p.root}, {p.energy()}")
            waiting_peaks = waiting_peaks[max_jump + 1:]
