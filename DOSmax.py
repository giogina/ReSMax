import re

from scipy import signal
import argparse
import numpy as np
from scipy.signal import argrelextrema

import concurrent.futures
import threading
import os
import subprocess

from DOSpeak import DOSpeak
from data_parser import parse, project_directory
import plot
from plot import resonance_summary_grid
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
        if type(root) is int:
            den = np.abs(data[root][2:] - data[root][:-2])
            data[f"rho_{root}"] = np.abs(num / den) if 0 not in list(den) else np.zeros(len(num))
    return data


def fitDOS(data, energy_range, thresholds, project_dir):
    """
    Fit DOS data using Lorentzian models and find resonances.

    Parameters:
    data (dict): Parsed data containing energy and DOS arrays.
    energy_range (tuple): Tuple specifying the energy range to consider for fitting.
    thresholds (list): List of energy thresholds for categorizing resonances.
    project_dir (str): Path to the project directory
    """
    result_file = project_dir + "resonances.txt"
    data_keys = list(data.keys())
    data_keys = [k for k in data_keys if type(k) is str and k.startswith("rho_")]
    data_keys.sort(key=lambda k: min(data[int(k[4:])]))  # It seems DOS maxima are best sampled as they come into view bottom->top.
    fitted_peaks_by_root = {}
    fitted_half_peaks_by_root = {}
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
        fitted_half_peaks_by_root[root] = []
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
            min_i = 0 if len(local_mins[0]) == 0 else max(local_mins[0])+3
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
                        dp.fit_lorentzian()  # todo: use the half slope idea to make better initial guesses
                        if dp.fit_E is not None and dp.energy() > lowest_populated_threshold and dp.warning is None:
                            fitted_peaks_by_root[root].append(dp)

    find_resonances(fitted_peaks_by_root)
    for res in Resonance.resonances:
        res.categorize_by_thresholds(thresholds)
    Resonance.resonances.sort(key=lambda r: r.energy)
                
                
def print_result_file(result_file):
    current_threshold = None
    with open(result_file, 'a') as save_file:
        for res in Resonance.resonances:
            if res.threshold != current_threshold:
                current_threshold = res.threshold
                print(f"\nResonances found below threshold E = {current_threshold}:")
                save_file.write(f"\r\nResonances found below threshold {current_threshold}:\r\n")
                save_file.write(f"Energy             Root\tgamma              \tSSR                \tRel. SSR per point\tGamma              \tA                 \ty0                \tOther roots\r\n")
            if res.best_fit is not None:
                print(f"{res.energy}, Root {res.best_fit.root}, SSR {res.best_fit.ssr}, Rel. SSR per point {res.best_fit.rel_ssr_per_point}, gamma = {res.best_fit.fit_gamma}, Gamma = {res.best_fit.fit_Gamma}, A = {res.best_fit.fit_A}, y0 = {res.best_fit.fit_y0}, Other roots = {[p.root for p in res.peaks]}")
                save_file.write(
                    f"{res.energy:.15f}".ljust(18, '0') + "\t" +
                    f"{res.best_fit.root}\t" +
                    f"{res.best_fit.fit_gamma:.15f}".ljust(18, '0') + "\t" +
                    f"{res.best_fit.ssr:.15f}".ljust(18, '0')[:18] + "\t" +
                    f"{res.best_fit.rel_ssr_per_point:.15f}".ljust(18, '0') + "\t" +
                    f"{res.best_fit.fit_Gamma:.15f}".ljust(18, '0') + "\t" +
                    f"{res.best_fit.fit_A:.15f}".ljust(18, '0') + "\t" +
                    f"{res.best_fit.fit_y0:.15f}".ljust(18, ' ') + "\t" +
                    f"{[p.root for p in res.peaks]} \t" +
                    f"{'[!] '+res.best_fit.warning if res.best_fit.warning is not None else ''}"
                    +"\r\n"
                )

        save_file.write("\r\n\r\n\r\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\r\n\r\n  Detailed data on all DOS peak fits:  \r\n\r\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\r\n\r\n")
        save_file.write(f"Threshold     \t\t Resonance         \t Energy (pointwise)\t Energy (fit)\t   Root\tgamma              \tSSR                \tRel. SSR per point\tGamma              \tA                 \ty0\r\n")

        current_threshold = None
        for res in Resonance.resonances:
            if res.threshold != current_threshold:
                current_threshold = res.threshold
            for peak in sorted(res.peaks, key=lambda p1: p1.root):
                save_file.write( " " +
                    f"{current_threshold:.3f}".ljust(10, ' ') + "\t" + (" x" if res.best_fit is not None and res.best_fit.root == peak.root else "  ") + "    \t" +
                    f"{res.energy:.15f}".ljust(18, '0') + "\t" +
                    f"{peak.pointwise_energy:.15f}".ljust(18, '0') + "\t" +
                    f"{peak.fit_E:.15f}".ljust(18, '0') + "\t" +
                    f"{peak.root}\t" +
                    f"{peak.fit_gamma:.15f}".ljust(18, '0') + "\t" +
                    f"{peak.ssr:.15f}".ljust(18, '0')[:18] + "\t" +
                    f"{peak.rel_ssr_per_point:.15f}".ljust(18, '0') + "\t" +
                    f"{peak.fit_Gamma:.15f}".ljust(18, '0') + "\t" +
                    f"{peak.fit_A:.15f}".ljust(18, '0') + "\t" +
                    f"{peak.fit_y0:.15f}".ljust(21, ' ') + "\t" +
                    f"{'[!] ' + peak.warning if peak.warning is not None else ''}"
                    + "\r\n"
                )

    print(f"\nResults have been written to {result_file}.")


def sliding_window_average(data, e_min, e_max, window_size):
    """
    Compute an averaged curve from a sliding window of energy-consecutive points using numpy.

    Parameters:
    data (dict): Parsed data dictionary containing "gamma", "rho_{root}" arrays, and others.
    e_min (float): Minimum energy value to include in the analysis.
    e_max (float): Maximum energy value to include in the analysis.
    window_size (int): Number of consecutive points in each sliding window.

    Returns:
    tuple: (averaged_energy, averaged_dos) where:
        - averaged_energy is the array of averaged energy values.
        - averaged_dos is the corresponding array of averaged DOS values.
    """
    # Collect energy and DOS points from the data
    points = []
    for key in data.keys():
        if isinstance(key, str) and key.startswith("rho_"):
            rho = data[key]
            energy = data[int(key[4:])]
            points.extend(zip(energy[1:-1], np.log10(rho)))  # Skip endpoints

    points = np.array(points)

    # Filter points within the energy range
    points = points[(points[:, 0] >= e_min) & (points[:, 0] <= e_max)]

    # Sort points by increasing energy
    points = points[np.argsort(points[:, 0])]

    # Apply a sliding window using numpy's stride tricks
    energy = points[:, 0]
    dos = points[:, 1]

    shape = (energy.size - window_size + 1, window_size)
    strides = (energy.strides[0], energy.strides[0])

    energy_windows = np.lib.stride_tricks.as_strided(energy, shape=shape, strides=strides)
    dos_windows = np.lib.stride_tricks.as_strided(dos, shape=shape, strides=strides)

    # Compute averages for each window
    averaged_energy = energy_windows.mean(axis=1)
    averaged_dos = dos_windows.mean(axis=1)

    return averaged_energy, averaged_dos, []


def filter_peaks_by_relative_prominence(clustering_array, peak_indices, prominence_threshold=0.3):
    """
    Filter peaks by relative prominence, processing from lowest to highest height.

    Parameters:
    clustering_array (array): The array of clustering weights.
    peak_indices (array): Indices of detected peaks in the clustering array.
    prominence_threshold (float): Minimum relative prominence required for a peak.

    Returns:
    array: Filtered peak indices that pass the relative prominence check.
    """
    # Sort peaks by their height in ascending order
    peak_heights = clustering_array[peak_indices]
    sorted_peak_indices = [idx for _, idx in sorted(zip(peak_heights, peak_indices))]

    # Convert to a set for efficient removal
    valid_peaks = set(sorted_peak_indices)

    for current_peak in sorted_peak_indices:
        if current_peak not in valid_peaks:
            continue  # Skip peaks already removed

        left_min = None
        for i in range(current_peak - 1, -1, -1):  # Search left
            if i in valid_peaks and i < current_peak:
                break  # Stop at the nearest valid left peak
            if left_min is None or clustering_array[i] < left_min:
                left_min = clustering_array[i]

        right_min = None
        for i in range(current_peak + 1, len(clustering_array)):  # Search right
            if i in valid_peaks and i > current_peak:
                break  # Stop at the nearest valid right peak
            if right_min is None or clustering_array[i] < right_min:
                right_min = clustering_array[i]

        # Calculate the prominence relative to the larger adjacent minimum
        if left_min is not None and right_min is not None:
            min_adjacent = max(left_min, right_min)
        elif left_min is not None:
            min_adjacent = left_min
        elif right_min is not None:
            min_adjacent = right_min
        else:
            min_adjacent = 0  # Shouldn't happen; safeguard

        relative_prominence = (clustering_array[current_peak] - min_adjacent) / clustering_array[current_peak]

        if relative_prominence < prominence_threshold:
            valid_peaks.remove(current_peak)

    return np.array(sorted(valid_peaks))


def cleanup(data):
    """
    Clean up the front of the data for each root by removing points before the first local minimum
    in the energy vs. gamma relationship.

    Parameters:
    data (dict): Parsed data dictionary containing "gamma", "{rootnumber}" arrays, and others.

    Returns:
    dict: Updated data dictionary with new keys "cleaned_{rootnumber}" containing valid indices.
    """

    cleaned_data = data.copy()
    for key in data.keys():
        if isinstance(key, int):
            energy = cleaned_data[key]
            cut_index = last_local_minimum(energy)
            valid_indices = list(range(cut_index, len(energy)))
            cleaned_data[f"cleaned_{key}"] = np.array(valid_indices)
    return cleaned_data


def last_local_minimum(energy):
    local_minima = argrelextrema(energy, np.less)[0]
    if len(local_minima) > 0:
        highest_local_min_index = local_minima[-1]
        local_min_energy = energy[highest_local_min_index]
        li, prev_lower_energy = next(
            ((i, e) for i, e in reversed(list(enumerate(energy[:highest_local_min_index]))) if e < local_min_energy),
            (None, None)  # Default value if no match is found
        )

        if prev_lower_energy is not None:
            hill_height = max(energy[li:highest_local_min_index])-local_min_energy
            if hill_height < 5*(10)**(-5) and (highest_local_min_index - li) < 4 and (local_min_energy - prev_lower_energy) < 10**(-4) and li not in local_minima:
                energy[li+1:highest_local_min_index] = np.linspace(energy[li], energy[highest_local_min_index], highest_local_min_index - li + 1)[1:-1]
                return last_local_minimum(energy)
        return highest_local_min_index + 3
    return 0


def main(file):
    """
    Main function to manage the DOS fitting process.

    Parameters:
    file (str): The path to the file to process.
    """
    print(f"Processing file: {file}")
    data = parse(file)
    # data = cleanup(data)  # todo: temp: no cleaning
    action = "o"  # overview plot
    low = None
    high = 0
    thresholds = None
    min_e = None
    max_e = None
    redo_overview = True
    while action != "r":
        if redo_overview:
            plot_file = project_directory(file) + "overview.png"
            min_e, max_e = plot.overview(data, plot_file, low, high)
            data = computeDOS(data)
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
                            "    'p': Plot panorama log(DOS) vs E\n"
                            # f"    'd': in case of discontinuity: Use {'point-wise maximum' if DOSpeak.discontinuity_treatment == 'fit' else 'fit'} (currently: using {DOSpeak.discontinuity_treatment})\n"
                            "    'x': exit\n"
                            )
        try:
            if next_action == 'rspp':
                criterion = 'rspp'
                print("DOS fits sorted by least relative SSR per data point.")
            elif next_action == 'ssr':
                criterion = 'ssr'
                print("DOS fits sorted by least sum of square residues (SSR).")
            # elif next_action == 'd':
            #     DOSpeak.discontinuity_treatment = 'point-wise maximum' if DOSpeak.discontinuity_treatment == 'fit' else 'fit'
            elif next_action == 'x':
                exit()
            elif next_action == 'p':
                plot.plot_all_resonance_peaks(data, Resonance.resonances, project_directory(file)[:-1] + "_dos_panorama.png")
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
                        thresholds.append(0)
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
                        if len(thresholds) and thresholds[-1] < 0:
                            thresholds.append(0)
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

    
    fitDOS(data, (low, high), thresholds, project_directory(file))

    # current_threshold = thresholds[0]
    # for res in Resonance.resonances:
    #     if res.threshold != current_threshold:
    #         plot.plot_resonance_partitions_with_clustering(data, Resonance.resonances, max(res.energy-5*res.best_fit.fit_Gamma, current_threshold), res.threshold, f"{plot.threshold_dir(project_directory(file), res.threshold)}all_res.png")
    #         current_threshold = res.threshold

    user_done = False
    max_thr = max([r.threshold for r in Resonance.resonances if r.threshold is not None])
    i = 1
    while i<len(thresholds):
        threshold = thresholds[i]
        resonance_overview_range = [thresholds[i-1], threshold]
        overview_plot_name = f"{plot.threshold_dir(project_directory(file), threshold)}all_resonances_{threshold:.3f}.png"
        if user_done:
            break

        if i > 0 and threshold <= max_thr:
            print(f"Plotting resonance overview for threshold {threshold}...")
            plot.plot_resonance_partitions_with_clustering(data, Resonance.resonances, resonance_overview_range[0], resonance_overview_range[1], overview_plot_name)

            while True:
                action = input(f"\nPlease verify the detected resonances in {project_directory(file)}resonance_plots.\n"
                          f"    'ok': Accept & proceed to next threshold\n"
                        + (f"    'back': Return to previous threshold\n" if (i>1) else "") +
                          f"    'end': Skip remaining thresholds & print results\n"
                          f"    'iRj': For resonance #i, select the peak of root #j\n"
                          f"    'iR': De-select resonance #i\n"
                          f"    'grid i': Plot all DOS peaks for resonance #i\n"
                          f"    'plot Emin Emax': Create resonance overview plot for E=Emin..Emax\n").strip()
                if action.lower() == "ok": # todo: close plots opened during this loop
                    i = i + 1
                    break
                elif action.lower() == "back":
                    if i > 1:
                        i = i - 1
                    break
                elif action.lower() == "end": # todo: do not print resonances after end
                    user_done = True
                    break
                elif action.lower().startswith("grid"): # todo: mark active plot in grid
                    _, i = action.split()
                    print("grid", i)
                    resonance_summary_grid(project_directory(file), Resonance.resonances, int(i))
                elif action.startswith("plot"):
                    try:
                        _, emin, emax = action.split()
                        resonance_overview_range = [float(emin), float(emax)]
                        overview_plot_name = f"{plot.threshold_dir(project_directory(file), threshold)}resonances_{emin}_{emax}.png"
                        print("plot", resonance_overview_range)

                        plot.plot_resonance_partitions_with_clustering(data, Resonance.resonances, resonance_overview_range[0], resonance_overview_range[1], overview_plot_name)
                    except ValueError:
                        print("Invalid format. Use: plot Emin Emax, e.g. plot -0.7 -0.5")
                else:
                    changes = re.findall(r'(\d+)R(\d+)?', action)
                    if not len(changes):
                        print("Invalid input.")
                    changed_thresholds = []
                    for change in changes:
                        res_index = int(change[0])
                        if res_index < len(Resonance.resonances):
                            res = Resonance.resonances[res_index]
                            if not res.threshold == threshold:
                                print(f"Resonance {res.index} at E={res.energy:.5f} does not belong to current threshold {threshold}, and is therefore skipped.")
                            else:
                                if change[1]:
                                    root_peaks = [p for p in res.peaks if p.root == int(change[1])]
                                    if len(root_peaks):
                                        res.best_fit = root_peaks[0]
                                        res.energy = root_peaks[0].fit_E
                                        changed_thresholds.append(res.threshold)
                                    else:
                                        print(f"Root {change[1]} does not contribute to Resonance {change[0]} at E={res.energy:.5f}; the change {change[0]}R{change[1]} is therefore skipped.")
                                else:
                                    res.best_fit = None
                                    changed_thresholds.append(res.threshold)
                    if i > 0 and threshold in changed_thresholds: # redo overview plot
                        plot.plot_resonance_partitions_with_clustering(data, Resonance.resonances, resonance_overview_range[0], resonance_overview_range[1], overview_plot_name)

# TODO: adjust spacing of 3R33 annotations; make active one bold.

    print_result_file(project_directory(file) + "resonances.txt")
    #  Todo:
    #   * implement toghether-fitting of double peaks
    #   * Better estimate of initial fit parameters

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
