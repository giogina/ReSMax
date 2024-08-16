from DOSpeak import DOSpeak, verbose


class Resonance:

    resonances = []

    def __init__(self, peak):
        """
        Initialize the Resonance object.

        Parameters:
        peak (DOSpeak or list): A DOSpeak object or a list of DOSpeak objects associated with this resonance.
        """
        if type(peak) is list:
            self.peaks = peak
            fit_qualities = [p.rel_ssr_per_point for p in peak]
            self.best_fit = self.peaks[fit_qualities.index(min(fit_qualities))]
            self.energy = self.best_fit.energy()
        else:
            self.peaks = [peak]
            self.energy = peak.energy()
            self.best_fit = peak
        self.resonances.append(self)
        self.threshold = None

    def add_peak(self, peak: DOSpeak):
        """
        Add a DOSpeak to the resonance and update the best fit if necessary.

        Parameters:
        peak (DOSpeak): The DOSpeak object to add.
        """
        self.peaks.append(peak)
        if peak.rel_ssr_per_point < self.best_fit.rel_ssr_per_point:
            self.best_fit = peak
            self.energy = peak.energy()

    def categorize_by_thresholds(self, thresholds):
        """
        Categorize the resonance by given energy thresholds.

        Parameters:
        thresholds (list): List of threshold energy values.
        """
        above = [t for t in thresholds if t > self.energy]
        if len(above):
            self.threshold = min(above)


def closest_resonance(resonances, peak: DOSpeak):
    """
    Find the closest resonance to a given DOSpeak.

    Parameters:
    resonances (list): List of existing Resonance objects.
    peak (DOSpeak): The DOSpeak object to compare.

    Returns:
    tuple: The index of the closest resonance and the error value, or (None, None) if no close resonance is found.
    """
    min_error = abs(peak.energy())
    closest_res_i = None
    for i, res in enumerate(resonances):
        if peak.root not in [p.root for p in res.peaks] and abs(res.energy - peak.energy()) < min_error:
            min_error = abs(res.energy - peak.energy())
            closest_res_i = i
    if closest_res_i is None:
        return None, None
    else:
        return closest_res_i, min_error


def assign_root_peaks_to_resonances(resonances, peaks):
    """
    Recursively assign root peaks to the closest resonance.

    Parameters:
    resonances (list): List of existing Resonance objects.
    peaks (list): List of DOSpeak objects to assign.
    """
    min_rel_err = 100
    closest_res = None
    closest_peak = None
    for i, peak in enumerate(peaks):
        cl, err = closest_resonance(resonances, peak)
        if cl is not None and err < min_rel_err:  # select closest well-matched peak
            min_rel_err = err * peak.rel_ssr_per_point
            closest_res = cl
            closest_peak = i

    if closest_res is None:
        for peak in peaks:
            print(f"New resonance: Root #{peak.root} {peak.energy()}, RSSRPP = {peak.rel_ssr_per_point}")
            Resonance(peak)  # new resonance energies found
    else:
        print(f"Assigning: Resonance {resonances[closest_res].energy}, peak {peaks[closest_peak].energy()} of root {peaks[closest_peak].root}, RSSRPP = {peaks[closest_peak].rel_ssr_per_point}")
        resonances[closest_res].add_peak(peaks[closest_peak])
        assign_root_peaks_to_resonances(resonances[:closest_res], peaks[:closest_peak])
        assign_root_peaks_to_resonances(resonances[closest_res + 1:], peaks[closest_peak + 1:])


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
            dists = [waiting_peaks[i + 1].energy() - waiting_peaks[i].energy() for i in range(0, len(waiting_peaks) - 1)]
            max_jump = dists.index(max(dists))
            res = Resonance(waiting_peaks[:max_jump + 1])
            if verbose:
                print("~~~~~", res.energy)
                for p in waiting_peaks[:max_jump + 1]:
                    print(f"    {p.root}, {p.energy()}")
            waiting_peaks = waiting_peaks[max_jump + 1:]

