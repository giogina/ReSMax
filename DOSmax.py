import os
from scipy import signal
import argparse
import numpy as np

from DOSpeak import DOSpeak
from data_parser import parse, project_directory
import plot
from resonance import find_resonances, Resonance

verbose = False

def computeDOS(data):
    """
    Compute the Density of States (DOS) based on gamma and root data.

    Parameters:
    data (dict): Parsed data containing gamma and root arrays.

    Returns:
    dict: Data with additional DOS arrays added.
    """
    roots = list(data.keys())
    roots.remove("gamma")
    num = data["gamma"][2:] - data["gamma"][:-2]  # This seems unnecessarily inaccurate. Why not i to i+1, with E in the middle?
    for root in roots:
        den = (data[root][2:] - data[root][:-2])
        data[f"rho_{root}"] = np.abs(num / den) if 0 not in list(den) else np.zeros(len(num))
    return data


def fitDOS(data, energy_range, fit_criterion, thresholds, result_file):
    """
    Fit DOS data using Lorentzian models and find resonances.

    Parameters:
    data (dict): Parsed data containing energy and DOS arrays.
    energy_range (tuple): Tuple specifying the energy range to consider for fitting.
    fit_criterion (str): Criterion for fitting ('ssr' or 'rspp').
    thresholds (list): List of energy thresholds for categorizing resonances.
    result_file (str): Path to the file where results will be written.
    """
    data_keys = list(data.keys())
    data_keys = [k for k in data_keys if type(k) is str and k.startswith("rho_")]
    data_keys.sort(
        key=lambda k: min(data[int(k[4:])]))  # It seems DOS maxima are best sampled as they come into view bottom->top.
    fitted_peaks_by_root = {}
    lowest_populated_threshold = None
    with open(result_file, 'w') as save_file:
        save_file.write("MBS:\r\n")
    print("Scanning all local DOS maxima...")
    for key in data_keys:
        root = int(key[4:])
        energy_array = data[root][1:-1]
        gamma_array = data["gamma"][1:-1]
        dos_array = data[key]
        fitted_peaks_by_root[root] = []
        thresholds_above = [t for t in thresholds if min(energy_array) < t]
        if len(thresholds_above):
            if lowest_populated_threshold is None:
                lowest_populated_threshold = thresholds_above[0]
            if lowest_populated_threshold == thresholds_above[0]:
                min_index = np.argmin(energy_array)
                mbs_string = f"Root {root}, E = {min(energy_array)} at gamma = {gamma_array[np.argmin(energy_array)]:.4f}"
                if min_index < 5 or min_index > len(gamma_array)-6:
                    mbs_string += f" [!] Minimum near {'lower' if min_index<5 else 'upper'} end of gamma range!"
                print(f"MBS: {mbs_string}")
                with open(result_file, 'a') as save_file:
                    save_file.write(f"{mbs_string}\r\n")
        if energy_range != (None, None):
            min_E = min(energy_array) if energy_range[0] is None else energy_range[0]
            max_E = max(energy_array) if energy_range[1] is None else energy_range[1]
            local_mins = signal.find_peaks(-energy_array)
            min_i = 0 if len(local_mins[0]) == 0 else max(local_mins[0])
            included_indices = [i for i, e in enumerate(energy_array) if (min_E <= e <= max_E) and min_i <= i]
            energy_array = np.array([energy_array[i] for i in included_indices])
            gamma_array = np.array([gamma_array[i] for i in included_indices])
            dos_array = np.array([dos_array[i] for i in included_indices])

        peaks, _ = signal.find_peaks(dos_array)
        valleys, _ = signal.find_peaks(-dos_array)
        if len(peaks) == 0:
            continue
        if len(valleys) == 0 or valleys[0] > peaks[0]:
            valleys = np.concatenate([np.array([0]), valleys])
        if len(valleys) == 0 or valleys[-1] < peaks[-1]:
            valleys = np.concatenate([valleys, np.array([len(dos_array) - 1])])

        if len(valleys) and len(peaks):
            ip = 0
            for iv in range(0, len(valleys) - 1):
                v = valleys[iv]
                v2 = valleys[iv + 1] + 1
                p = peaks[ip]
                while p < v and len(peaks) > ip + 1:
                    ip += 1
                    p = peaks[ip]
                if v < p < v2:  # successfully identified section
                    peak_E = float(energy_array[p])
                    peak_rho = float(dos_array[p])
                    if v2 - v >= 10:  # necessary to fit all parameters
                        energies = energy_array[v:v2]
                        gammas = gamma_array[v:v2]
                        rhos = dos_array[v:v2]
                        dp = DOSpeak(energies, rhos, gammas, root, peak_E, peak_rho)
                        fitted_paras = dp.fit_lorentzian()
                        if fitted_paras is not None and dp.energy() > lowest_populated_threshold:
                            fitted_peaks_by_root[root].append(dp)

    find_resonances(fitted_peaks_by_root)
    for res in Resonance.resonances:
        res.categorize_by_thresholds(thresholds)
    Resonance.resonances.sort(key=lambda res: res.energy)

    current_threshold = None
    with open(result_file, 'a') as save_file:
        for res in Resonance.resonances:
            if res.threshold != current_threshold:
                current_threshold = res.threshold
                print(f"\nResonances found below threshold {current_threshold}:")
                save_file.write(f"\r\nResonances found below threshold {current_threshold}:\r\n")
                save_file.write(f"Energy             Root\tgamma              \tSSR                \tRel. SSR per point\tGamma              \tA                 \ty0                \tOther roots\r\n")
            print(f"{res.energy}, Root {res.best_fit.root}, SSR {res.best_fit.ssr}, Rel. SSR per point {res.best_fit.rel_ssr_per_point}, gamma = {res.best_fit.fit_gamma}, Gamma = {res.best_fit.fit_Gamma}, A = {res.best_fit.fit_A}, y0 = {res.best_fit.fit_y0}, Other roots = {[p.root for p in res.peaks]}")
            save_file.write(
                f"{res.energy:.15f}".ljust(18, '0') + "\t" +
                f"{res.best_fit.root}\t" +
                f"{res.best_fit.fit_gamma:.15f}".ljust(18, '0') + "\t" +
                f"{res.best_fit.ssr:.15f}".ljust(18, '0')[:18] + "\t" +
                f"{res.best_fit.rel_ssr_per_point:.15f}".ljust(18, '0') + "\t" +
                f"{res.best_fit.fit_Gamma:.15f}".ljust(18, '0') + "\t" +
                f"{res.best_fit.fit_A:.15f}".ljust(18, '0') + "\t" +
                f"{res.best_fit.fit_y0:.15f}".ljust(18, '0') + "\t" +
                f"{[p.root for p in res.peaks]} \t" +
                f"{'!Warning! '+res.best_fit.warning if res.best_fit.warning is not None else ''}"
                +"\r\n"
            )
    print(f"\nResults have been written to {result_file}.")

    # TODO: better treatment of pointwise maximum energy - e.g. the figures are still named using fit_E I think? Also, show whole range (including dropped points) of fit.
    #  Also, why is first [19] not included in following resonance? Use other E for that?
    #  Allow turning on/off of auto-trim.
    #  Give a slightly larger allowance for the fit_check.


def main(file):
    """
    Main function to manage the DOS fitting process.

    Parameters:
    file (str): The path to the file to process.
    """
    print(f"Processing file: {file}")
    data = parse(file)
    action = "o"  # overview plot
    low = None
    high = 0
    criterion = 'rspp'
    thresholds = None
    min_e = None
    max_e = None
    redo_overview = True
    while action != "r":
        if redo_overview:
            plot_file = project_directory(file) + "overview.png"
            min_e, max_e = plot.overview(data, plot_file, low, high)
            if thresholds is None:
                thresholds = [max_e]  # dummy value in case no thresholds are entered
            redo_overview = False
            print(f"\nγ vs E overview graph has been plotted to {plot_file}.")
        next_action = input(f"Please specify your next action: \n\n"
                            "    'o': Plot gamma vs E overview graph (full energy range).\n"
                            "    'o E_min E_max': Re-plot gamma vs E overview graph (specified energy range).\n"
                            "    'r': Fit density of state over the currently displayed range ({'-∞' if low is None else low}..{high})\n"
                            "    'r E_min E_max': Fit density of state in the range E_min..E_max\n"
                            "    'z': input nuclear charge Z\n"
                            "    't': input list of thresholds\n"
                            f"    'd': in case of discontinuity: Use {'point-wise maximum' if DOSpeak.discontinuity_treatment == 'fit' else 'fit'} (currently: using {DOSpeak.discontinuity_treatment})\n"
                            "    'x': exit\n"
                            )
        try:
            if next_action == 'rspp':
                criterion = 'rspp'
                print("DOS fits sorted by least relative SSR per data point.")
            elif next_action == 'ssr':
                criterion = 'ssr'
                print("DOS fits sorted by least sum of square residues (SSR).")
            elif next_action == 'd':
                DOSpeak.discontinuity_treatment = 'point-wise maximum' if DOSpeak.discontinuity_treatment == 'fit' else 'fit'
            elif next_action == 'x':
                exit()
            elif next_action == 'z':
                ok = False
                while not ok:
                    inp_z = input("Z = ")
                    try:
                        z = int(inp_z)
                        n = 1
                        thresholds = [-z * z / 2.]
                        if max_e is None or max_e >= 0:
                            max_e = -z * z / 2. / 400  # should never happen - limit to 20 threshold values
                        while thresholds[-1] < max_e:
                            n += 1
                            thresholds.append(-z * z / 2. / n / n)
                        print(f"Threshold values: {' '.join([str(t) for t in thresholds if t > min_e])}")
                        ok = True
                    except:
                        print("Please input an integer.")
            elif next_action == 't':
                ok = False
                while not ok:
                    inp_t = input("Threshold values (separate by space) = ")
                    try:
                        thresholds = [float(t) for t in inp_t.split()]
                        thresholds.sort()
                        print(f"Threshold values: {' '.join([str(t) for t in thresholds])}")
                        ok = True
                    except Exception as e:
                        print("Please input a space-separated list of floats.")
            else:
                response = next_action.split(' ')
                if len(response) == 3:
                    action = response[0]
                    low = float(response[1])
                    high = float(response[2])
                    if action == "o":
                        redo_overview = True
                elif len(response) == 1 and response[0] == "o":  # reset overview
                    low = None
                    high = 0
                    redo_overview = True
                elif len(response) == 1 and response[0] == "r":  # proceed to DOS fit
                    break
                else:
                    print("Invalid input format. Please input action in the form e.g. 'o -0.2 0'.")
        except Exception as e:
            print(e)
            print("Invalid input format. Please input action in the form e.g. 'o -0.2 0'.")

    data = computeDOS(data)
    fitDOS(data, (low, high), criterion, thresholds, result_file=project_directory(file) + "resonances.txt")

    input_string = f"\nPlease specify your next action: \n\n    'p': plot best fit for all resonances\n"
    populated_thresholds = set([res.threshold for res in Resonance.resonances])
    for i, t in enumerate(thresholds):
        if t in populated_thresholds:
            input_string += f"    'p{i}': plot best fit for all resonances below threshold {t}\n"
    input_string += "    'x': exit\n"
    while True:
        action = input(input_string)
        if action == 'p':
            plot.resonance_fits(project_directory(file), Resonance.resonances)
            break
        elif action in [f'p{i}' for i, t in enumerate(thresholds)]:
            i = int(action[1:])
            plot.resonance_fits(project_directory(file), Resonance.resonances, thresholds[i])
        elif action == 'x':
            break
        else:
            print("Invalid input.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a file.")
    parser.add_argument('-f', '--file', type=str, required=True, help="Path to the file")
    args = parser.parse_args()
    main(args.file)
