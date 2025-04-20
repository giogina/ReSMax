import cProfile
import pstats
import sys
import os
import subprocess
import importlib.util
import importlib.metadata
import time

verbose = False

# Freeze: pyinstaller --onefile DOSmax.py

def install_requirements():
    """Ensure all required dependencies are installed."""

    # Check if running as a frozen .exe
    if getattr(sys, 'frozen', False):
        return

    # Check for requirements.txt
    req_file = "requirements.txt"
    if not os.path.exists(req_file):
        print("No requirements.txt found. Skipping dependency check.")
        return

    # Read required packages (ignore version specifiers for easier installation)
    with open(req_file, "r") as f:
        required = [line.strip().split("==")[0] for line in f if line.strip()]

    # Get installed packages
    installed = {pkg.metadata["Name"].lower() for pkg in importlib.metadata.distributions()}
    missing = [pkg for pkg in required if pkg.lower() not in installed]

    # Install missing dependencies
    if missing:
        print(f"Installing missing dependencies: {', '.join(missing)}, please wait...")

        try:
            # Skip upgrade if pip is at least 22.0
            current_version = importlib.metadata.version("pip")
            if tuple(map(int, current_version.split("."))) < (22, 0):
                print("Upgrading pip...")
                subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "pip"])
            import pip

        except ImportError:
            print("pip is not installed. Attempting to install pip...")
            try:
                subprocess.check_call([sys.executable, "-m", "ensurepip"])
                subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "pip"])
                import pip
            except Exception as e:
                print(f"Failed to install pip: {e}")
                return

        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", *missing])
            print("All required packages installed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error installing packages: {e}")

install_requirements()

from scipy import signal
import argparse
import numpy as np
from scipy.signal import argrelextrema
import re

from DOSpeak import DOSpeak
from data_parser import parse, project_directory
import plot
from resonance import find_resonances, Resonance


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
            den = data[root][2:] - data[root][:-2]
            data[f"rho_{root}"] = (num / den) if 0 not in list(den) else np.zeros(len(num))
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
    print("Scanning all local DOS maxima...\n")
    valley_points = []
    save_file_string = "MBS:\r\n"

    min_energy = np.min([np.min(data[root][1:-1]) for root in data if isinstance(root, int)])
    max_energy = np.max([np.max(data[root][1:-1]) for root in data if isinstance(root, int)])

    min_E = min_energy if energy_range[0] is None else energy_range[0]  # for masking data
    max_E = max_energy if energy_range[1] is None else energy_range[1]

    lowest_populated_threshold = ([t for t in thresholds if min_energy < t] + [0])[0]

    for key in data_keys:
        root = int(key[4:])
        energy_array = data[root][1:-1]
        gamma_array = data["gamma"][1:-1]
        min_e = np.min(energy_array)
        if min_e < lowest_populated_threshold:
            min_index = np.argmin(energy_array)
            mbs_string = f"Root {root}, E = {min_e} at gamma = {gamma_array[min_index]:.4f}"
            if min_index < 5 or min_index > len(gamma_array)-6:
                mbs_string += f" [!] Minimum near {'lower' if min_index<5 else 'upper'} end of gamma range!"
            print(f"MBS: {mbs_string}")
            save_file_string += f"{mbs_string}\r\n"
    with open(result_file, 'w') as save_file:
        save_file.write(save_file_string)

    print("\n Detecting and fitting DOS peaks...")

    for key in data_keys:
        root = int(key[4:])
        energy_array = data[root][1:-1]
        mask = (energy_array >= min_E) & (energy_array <= max_E)
        if not np.any(mask):
            continue
        energy_array = energy_array[mask]
        gamma_array = data["gamma"][1:-1][mask]
        dos_array = data[key][mask]
        if len(energy_array) < 10:
            continue
        deriv_array = np.gradient(energy_array, gamma_array)  # to prevent the two-index jump in DOS from obscuring things


        fitted_peaks_by_root[root] = []
        # valley_indices, _ = signal.find_peaks(-dos_array)
        valley_indices, _ = signal.find_peaks(deriv_array)

        valley_indices = [i for i in valley_indices if dos_array[i] >= 0]  # only cut on ascends, not in the middle of descending pieces.
        if len(valley_indices) < 2 or dos_array[1] > dos_array[0]:
            valley_indices = np.concatenate([np.array([0]), valley_indices])
        if len(valley_indices) < 2 or dos_array[-2] > dos_array[-1]:
            valley_indices = np.concatenate([valley_indices, np.array([len(dos_array) - 1])])

        for iv in range(0, len(valley_indices) - 1):
            v = int(valley_indices[iv])+1
            v2 = int(valley_indices[iv+1])
            energies = energy_array[v:v2]
            gammas = gamma_array[v:v2]
            rhos = dos_array[v:v2]
            # valley_points.append([float(gamma_array[v]), float(energy_array[v])])

            if len(rhos) < 10:  # too few points
                continue

            dp = DOSpeak(energies, rhos, gammas, root)  # descending sections can be added; for these, rel_ssr_per_point = 10**6, energy() = min(energies)
            dp.fit_lorentzian()
            if dp.energy() is not None and dp.energy() > lowest_populated_threshold and dp.warning is None:
                fitted_peaks_by_root[root].append(dp)

    # plot.plot_partitions(data, fitted_peaks_by_root, project_dir+"partitions.png", valley_points)  # also uncomment the valley_points.append above to activate red split points

    print(" Finding resonances...")
    find_resonances(fitted_peaks_by_root)
    for res in Resonance.resonances:
        res.categorize_by_thresholds(thresholds)
    Resonance.resonances.sort(key=lambda r: r.energy)

                
def print_result_file(max_threshold, result_file):
    current_threshold = None
    with open(result_file, 'a') as save_file:
        for res in Resonance.resonances:
            show = res.should_be_shown()
            if show is False:
                continue
            if res.threshold != current_threshold:
                current_threshold = res.threshold
                if current_threshold > max_threshold:
                    break
                print(f"\nResonances found below threshold E = {current_threshold}:")
                save_file.write(f"\r\nResonances found below threshold {current_threshold}:\r\n")
                save_file.write(f"Energy       \t   Root\tgamma              \tSSR                \tRel. SSR per point\tGamma              \tA                 \ty0                \tOther roots\r\n")
            if show:
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
                    f"{'[!] '+res.best_fit.warning if res.best_fit.warning is not None else ''}" +
                    "\r\n"
                )
            elif show is None:  # descending section
                print(
                    f"{res.energy}, Root {res.best_fit.root}  [!] Descending section; DOS peak could not be fitted. [!]  Other roots = {[p.root for p in res.peaks]}")
                save_file.write(
                    f"{res.energy:.15f}".ljust(18, '0') + "\t" +
                    f"{res.best_fit.root}\t" +
                    f"{res.best_fit.fit_gamma:.15f}".ljust(18, '0') + "\t" +
                    f"\t[!] Descending section; DOS peak could not be fitted. [!]\t\t\t\t" +
                    f"{[p.root for p in res.peaks]} \t" +
                    "\r\n"
                )

        save_file.write("\r\n\r\n\r\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\r\n\r\n  Detailed data on all DOS peak fits:  \r\n\r\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\r\n\r\n")
        save_file.write(f"Threshold     \t\t Resonance         \t Energy (pointwise)\t Energy (fit)\t   Root\tgamma              \tSSR                \tRel. SSR per point\tGamma              \tA                 \ty0\r\n")

        current_threshold = None
        for res in Resonance.resonances:
            if res.threshold != current_threshold:
                current_threshold = res.threshold
            for peak in sorted(res.peaks, key=lambda p1: p1.root):
                if peak.is_descending:
                    save_file.write(" " +
                        f"{current_threshold:.3f}".ljust(10, ' ') + "\t" + (
                            " x" if res.best_fit is not None and res.best_fit.root == peak.root else "  ") + "    \t" +
                        f"{res.energy:.15f}".ljust(18, '0') + "\t" +
                        f"{peak.approx_peak_E:.15f}".ljust(18, '0') + "\t" +
                        f"       --- ".ljust(18, ' ') + "\t" +
                        f"{peak.root}\t" +
                        f"{peak.fit_gamma:.15f}".ljust(18, '0') + "\t" +
                        f"\t[!] Energy is descending, therefore there isn't a fitable DOS peak [!]\t" +
                        "\r\n"
                    )
                else:
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
                        f"{'[!] ' + peak.warning if peak.warning is not None else ''}" +
                        "\r\n"
                    )

    print(f"\nResults have been written to {result_file}.")
    plot.open_file(result_file)


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

def generate_thresholds(z: int):
    n = 1
    thresholds = [-z * z / 2.]
    # if max_E is None or max_E >= 0:
    #     max_E = -z * z / 2. / 400  # should never happen - limit to 20 threshold values
    # while thresholds[-1] < max_E:
    while n < 20:
        n += 1
        thresholds.append(-z * z / 2. / n / n)
    thresholds.append(0)
    return thresholds

def main(file):
    """
    Main function to manage the DOS fitting process.

    Parameters:
    file (str): The path to the file to process.   # TODO: for cat-stab example, MBS are under the 0.5 threshold. Make sure first plot starts above.
    """
    display_timing = False
    if display_timing:
        start = time.time()

    print(f"Processing file: {file}")
    data = parse(file)

    if display_timing:
        end = time.time()
        print("Parsing done ", end-start)
        print("Number of roots: ", len(data.keys()) - 1)
        print("Number of gamma values: ", len(data["gamma"]))

    action = "o"  # overview plot
    low = None
    high = 0
    thresholds = None
    redo_overview = True
    overview_margin = 0.03
    while action != "r":
        if redo_overview:
            plot_file = project_directory(file) + "stabilization_diagram.png"
            plot.overview(data, plot_file, low, high, margin=overview_margin)
            data = computeDOS(data)
            if thresholds is None:
                thresholds = generate_thresholds(2)  # dummy value in case no thresholds are entered
            redo_overview = False
            print(f"\nγ vs E overview graph has been plotted to {plot_file}.")
        next_action = input(f"Please specify your next action: \n\n"
                            "    'o': Plot stabilization diagram (full energy range).\n"
                            "    'o E_min E_max': Plot stabilization diagram for E = E_min .. E_max.\n"
                            "    'z': input nuclear charge Z (default: Z = 2)\n"
                            "    't': input list of thresholds (default: [-Z^2/n^2/2, n = 1..20])\n"
                            "    'p': Plot clustering panorama (log(DOS) vs E)\n"
                            "    'p E_min E_max': Plot clustering panorama (log(DOS) vs E) for E = E_min..Emax\n"
                            "    'r': Fit density of state over the currently displayed range ({'-∞' if low is None else low}..{high})\n"
                            "    'r E_min E_max': Fit density of state in the range E_min..E_max\n"
                            "    'x': exit\n"
                            )
        try:
            if next_action.lower() == 'rspp':
                criterion = 'rspp'
                print("DOS fits sorted by least relative SSR per data point.")
            elif next_action.lower() == 'ssr':
                criterion = 'ssr'
                print("DOS fits sorted by least sum of square residues (SSR).")
            # elif next_action.lower() == 'd':
            #     DOSpeak.discontinuity_treatment = 'point-wise maximum' if DOSpeak.discontinuity_treatment == 'fit' else 'fit'
            elif next_action.lower() == 'x':
                exit()
            elif next_action.lower() == 'p':
                plot.plot_all_resonance_peaks(data, Resonance.resonances, project_directory(file) + "dos_panorama.png")
            elif next_action.lower() == 'z':
                ok = False
                while not ok:
                    inp_z = input("Z = ")
                    try:
                        z = int(inp_z)
                        thresholds = generate_thresholds(z)
                        print(f"Threshold values: {' '.join([str(t) for t in thresholds])}")
                        ok = True
                    except:
                        print("Please input an integer.")
            elif next_action.lower() == 't':
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
                        overview_margin = 0
                    if action == "p":
                        plot.plot_all_resonance_peaks(data, Resonance.resonances, project_directory(file) + f"dos_panorama_{low}_{high}.png", low, high)
                elif len(response) == 1 and response[0] == "o":  # reset overview
                    low = None
                    high = 0
                    overview_margin = 0.03
                    redo_overview = True
                elif len(response) == 1 and response[0] == "r":  # proceed to DOS fit
                    break
                else:
                    print("Invalid input format. Please input action in the form e.g. 'o -0.2 0'.")
        except Exception as e:
            print(e)
            print("Invalid input format. Please input action in the form e.g. 'o -0.2 0'.")
    plot.start_clustering_background_preparation(data, thresholds)

    if display_timing:
        start = time.time()

    fitDOS(data, (low, high), thresholds, project_directory(file))

    if display_timing:
        end = time.time()
        print("Fitting done ", end-start)
        start = end

    max_thr = max([r.threshold for r in Resonance.resonances if r.threshold is not None], default=0) # governs check loop and resonances.txt output cutoff!
    i = 1
    while i<len(thresholds):
        threshold = thresholds[i]
        resonance_overview_range = [thresholds[i-1], threshold]
        overview_plot_name = f"{plot.threshold_dir(project_directory(file), threshold)}all_resonances_{threshold:.3f}.png"
        if threshold > max_thr:
            break

        if i > 0 and threshold <= max_thr:
            print(f" Plotting resonance overview for threshold {threshold}...")
            plot.resonance_partitions_with_clustering(data, Resonance.resonances, resonance_overview_range[0], resonance_overview_range[1], overview_plot_name, None, threshold_above=threshold)
            manual_range = False

            while True:
                action = input(f"\nPlease verify the detected resonances in {project_directory(file)}resonance_plots.\n\n"
                          f"    'ok': Accept & proceed to next threshold\n"
                        + (f"    'back': Return to previous threshold\n" if (i>1) else "") +
                          f"    'end': Skip remaining thresholds & print results\n"
                          f"    'iRj': For resonance #i, select the peak of root #j\n"
                          f"    'iR': Toggle resonance #i on/off\n"
                          f"    'i+R': De-select all resonances from #i up to the threshold {threshold}\n"
                          f"    'grid i': Plot all DOS peaks for resonance #i\n"
                          f"    'plot Emin Emax': Create resonance overview plot for E=Emin..Emax\n"
                          f"    'close': Close image viewer\n"
                               ).strip()
                try:
                    if action.lower() == "ok" or action.lower() == "next":
                        i = i + 1
                        break
                    elif action.lower() == "back" or action.lower().startswith("prev"):
                        if i > 1:
                            i = i - 1
                        break
                    elif action.lower() == "end":
                        max_thr = threshold
                        i = i + 1
                        break
                    elif action.lower().startswith("grid"):
                        _, n = action.split()
                        plot.resonance_summary_grid(project_directory(file), Resonance.resonances, int(n)-1, None)
                    elif action.startswith("plot"):
                        try:
                            _, emin, emax = action.split()
                            resonance_overview_range = [float(emin), float(emax)]
                            resonance_overview_range.sort()
                            overview_plot_name = f"{plot.threshold_dir(project_directory(file), threshold)}resonances_{emin}_{emax}.png"
                            manual_range = True
                            plot.resonance_partitions_with_clustering(data, Resonance.resonances, resonance_overview_range[0], resonance_overview_range[1], overview_plot_name, None, manual_range=manual_range, threshold_above=threshold)
                        except ValueError:
                            print("Invalid format. Use: plot Emin Emax, e.g. plot -0.7 -0.5")
                    elif action.lower() == "close":
                        plot.close_files(project_directory(file))
                    else:
                        changes = re.findall(r'(\d+)(\+?R)(\d+)?', action)
                        if not len(changes):
                            print("Invalid input.")
                        changed_thresholds = []
                        for change in changes:
                            res_index = int(change[0])-1
                            operator = change[1]
                            root_index = int(change[2]) if change[2] else None

                            if operator == 'R':
                                if res_index < len(Resonance.resonances):
                                    res = Resonance.resonances[res_index]
                                    if res.threshold == threshold:
                                        if root_index is not None:
                                            root_peaks = [p for p in res.peaks if p.root == root_index]
                                            if len(root_peaks):
                                                res.best_fit = root_peaks[0]
                                                res.manual_peak_selection = True
                                                res.energy = root_peaks[0].energy()
                                                changed_thresholds.append(res.threshold)
                                            else:
                                                print(f"Root {change[1]} does not contribute to Resonance {change[0]} at E={res.energy:.5f}; the change {change[0]}R{change[1]} is therefore skipped.")
                                        else:
                                            if (res.should_be_shown() != False):
                                                res.best_fit = None
                                            else:
                                                res.manual_peak_selection = True
                                                res.select_best_fit()
                                            changed_thresholds.append(res.threshold)
                                    else:
                                        print(f"Resonance {res.index+1} at E={res.energy:.5f} does not belong to current threshold {threshold}, and is therefore skipped.")

                            elif operator == '+R':
                                for idx in range(res_index, len(Resonance.resonances)):
                                    res = Resonance.resonances[idx]
                                    if res.threshold != threshold:
                                        break
                                    res.best_fit = None
                                    changed_thresholds.append(res.threshold)

                        if i > 0 and threshold in changed_thresholds:  # redo overview plot if necessary
                            plot.resonance_partitions_with_clustering(data, Resonance.resonances, resonance_overview_range[0], resonance_overview_range[1], overview_plot_name, None, manual_range=manual_range, threshold_above=threshold)
                        
                        # changes = re.findall(r'(\d+)R(\d+)?', action)
                        # if not len(changes):
                        #     print("Invalid input.")
                        # changed_thresholds = []
                        # for change in changes:
                        #     res_index = int(change[0])
                        #     if res_index < len(Resonance.resonances):
                        #         res = Resonance.resonances[res_index]
                        #         if not res.threshold == threshold:
                        #             print(f"Resonance {res.index} at E={res.energy:.5f} does not belong to current threshold {threshold}, and is therefore skipped.")
                        #         else:
                        #             if change[1]:
                        #                 root_peaks = [p for p in res.peaks if p.root == int(change[1])]
                        #                 if len(root_peaks):
                        #                     res.best_fit = root_peaks[0]
                        #                     res.energy = root_peaks[0].energy()
                        #                     changed_thresholds.append(res.threshold)
                        #                 else:
                        #                     print(f"Root {change[1]} does not contribute to Resonance {change[0]} at E={res.energy:.5f}; the change {change[0]}R{change[1]} is therefore skipped.")
                        #             else:
                        #                 res.best_fit = None
                        #                 changed_thresholds.append(res.threshold)
                        # if i > 0 and threshold in changed_thresholds: # redo overview plot
                        #     plot.resonance_partitions_with_clustering(data, Resonance.resonances, resonance_overview_range[0], resonance_overview_range[1], overview_plot_name, None, manual_range=manual_range, threshold_above=threshold)
                except Exception as e:
                    print(f"Invalid input: {e}")

    print_result_file(max_thr, project_directory(file) + "resonances.txt")
    #  Ideas:
    #   * Better estimate of initial fit parameters - use the half slope idea to make better initial guesses (see experimental_branch folder on laptop)

    input_string = f"\nPlease specify your next action: \n\n    'p': plot best fit for all resonances\n"
    populated_thresholds = set([res.threshold for res in Resonance.resonances])
    for i, t in enumerate(thresholds):
        if t in populated_thresholds:
            input_string += f"    'p{i}': plot best fit for all resonances below threshold {t}\n"
    input_string += "    'c': Close image viewer\n"
    input_string += "    'x': exit\n"
    while True:
        action = input(input_string)
        try:
            if action.lower() == 'p':
                plot.resonance_fits(project_directory(file), Resonance.resonances)
                break
            elif action.lower() in [f'p{i}' for i, t in enumerate(thresholds)]:
                i = int(action[1:])
                plot.resonance_fits(project_directory(file), Resonance.resonances, thresholds[i])
            elif action.lower().startswith('c'):
                plot.close_files(project_directory(file))
            elif action.lower() == 'x':
                break
            else:
                print("Invalid input.")
        except Exception as e:
            print(f"Invalid input: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a file.")
    parser.add_argument('-f', '--file', type=str, required=True, help="Path to the file")
    args = parser.parse_args()
    main(args.file)
    # with cProfile.Profile() as pr:
    #     main(args.file)
    # stats = pstats.Stats(pr)
    # stats.sort_stats("cumulative").print_stats(100)