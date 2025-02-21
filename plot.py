import os
import subprocess
import platform
import time

import numpy as np
import psutil

from concurrent.futures import ThreadPoolExecutor

from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, hsv_to_rgb
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import Affine2D



# matplotlib.use('Agg')  # faster backend

def get_root_color(rootnr: int):
    saturation = 1
    value = 0.8
    hue = (rootnr * 222.5+15) % 360
    hue /= 360.0
    return hsv_to_rgb((hue, saturation, value))


def plot_DOS(data, root, file, fitted_peaks_by_root=None):
    """
    Plot the Density of States (DOS) for a specific root and save the plot.

    Parameters:
    data (dict): The parsed data containing energy and DOS arrays.
    root (int): The root identifier for the data to plot.
    fitted_peaks_by_root (dict, optional): Dictionary of fitted peaks by root. Defaults to None.
    """
    x_data = data["gamma"][1:-1]
    y_data = data[root][1:-1]
    plt.figure(figsize=(12, 8))
    plt.scatter(x_data, y_data, edgecolor='black', facecolor='none')
    if fitted_peaks_by_root is not None:
        for peak in fitted_peaks_by_root[root]:
            plt.plot([min(x_data), max(x_data)], [peak.energy(), peak.energy()], 'r-')
    plt.savefig(file)
    plt.close()
    open_file(file)


def peak_fit(dos_peak, file):
    """
    Plot the fitted peak for a DOSpeak object and save the plot.

    Parameters:
    dos_peak (DOSpeak): The DOSpeak object containing the fitted data.
    file (str): The path where the plot should be saved.
    """
    x_data = dos_peak.energy_array
    y_data = dos_peak.dos_array
    xmin, xmax = dos_peak.energy()-8*dos_peak.fit_Gamma, dos_peak.energy()+8*dos_peak.fit_Gamma
    x_smooth = np.linspace(xmin, xmax, 1000)
    y_smooth = dos_peak.get_smooth_lorentzian_curve(x_smooth)
    plt.figure(figsize=(12, 8))
    plt.scatter(x_data, y_data, edgecolor="black", facecolor='none')
    plt.plot(x_smooth, y_smooth, 'r-')
    plt.xlabel("Energy (a.u.)")
    plt.ylabel("DOS")
    plt.xlim(xmin, xmax)
    plt.savefig(file)
    plt.close()
    with open(f"{file[:-5]}.txt", 'w') as plot_data_file:
        plot_data_file.write(np.array2string(x_data, separator=" ", max_line_width=np.inf) + "\r\n")
        plot_data_file.write(np.array2string(y_data, separator=" ", max_line_width=np.inf) + "\r\n")
        plot_data_file.write(np.array2string(x_smooth, separator=" ", max_line_width=np.inf) + "\r\n")
        plot_data_file.write(np.array2string(y_smooth, separator=" ", max_line_width=np.inf) + "\r\n")


def overview(data, plot_file, from_e=None, to_e=None, margin = 0.03):
    """
    Plot an overview of gamma vs energy and save the plot.

    Parameters:
    data (dict): The parsed data containing energy and gamma arrays.
    plot_file (str): The path where the plot should be saved.
    from_e (float, optional): Minimum energy value for the plot. Defaults to None.
    to_e (float, optional): Maximum energy value for the plot. Defaults to 0.

    Returns:
    tuple: Minimum and maximum energy values plotted.
    """
    plt.figure(figsize=(32, 24))
    min_E = 0 if from_e is None else from_e
    max_E = -100 if to_e is None else to_e

    for key, values in data.items():
        if type(key) is int:
            plt.plot(data["gamma"], values, label=key, color=get_root_color(key))  # plot / scatter
            if from_e is None:
                min_E = min(min_E, min(values))
            if to_e is None:
                max_E = max(max_E, max(values))

    plt.xlabel("γ")
    plt.ylabel("Energy (a.u.)")
    plt.ylim(bottom=min_E - margin, top=max_E + margin)
    plt.locator_params(axis='x', nbins=20)
    plt.locator_params(axis='y', nbins=20)
    plt.minorticks_on()

    plt.savefig(plot_file)
    plt.close()
    open_file(plot_file)

def open_file(file, opened_files = None):
    """
    Open a file using the default application based on the operating system.

    Parameters:
    file (str): The file path to open.
    """
    # (opened_files is not a useful list due to xdg-open automatically choosing a software. Could track pid's, but that's too complicated.)
    if platform.system() == 'Windows':
        os.startfile(file)
        # subprocess.Popen(["start", "", file], shell=True, creationflags=0x00000008, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif platform.system() == 'Darwin':  # macOS
        subprocess.call(('open', file), stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)  # .call necessary for .txt?
        # subprocess.Popen(["open", file], preexec_fn=os.setsid, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # proc = subprocess.Popen(['open', file])
    else:  # Linux and other Unix-like
        subprocess.call(('xdg-open', file), stderr=subprocess.DEVNULL)
        # subprocess.Popen(["xdg-open", file], preexec_fn=os.setsid, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # proc = subprocess.Popen(['xdg-open', file])
    # if opened_files is not None:
    #     opened_files.append(proc)  # Store process reference


def close_files(output_path):
    own_pids = []

    for proc in psutil.process_iter(["pid", "cmdline"]):
        try:
            cmdline = proc.info["cmdline"]
            if cmdline and any(output_path in arg for arg in cmdline):
                own_pids.append(proc.info["pid"])
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    if not len(own_pids):
        print("Image viewer process could not be identified. This only works if the viewer has been started by the current script.")
    for pid in own_pids:
        try:
            proc = psutil.Process(pid)
            proc.terminate()
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass


def threshold_dir(project_dir, threshold):
    sep = '\\' if '\\' in project_dir else '/'
    plots_dir = f"{project_dir}resonance_plots{sep}"
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    thres_dir = f"{plots_dir}{threshold:.5f}{sep}"
    if not os.path.exists(thres_dir):
        os.mkdir(thres_dir)
    return thres_dir


def resonance_fits(project_dir, resonances, threshold=None):
    """
    Plot the resonance fits and save them in the appropriate directories.

    Parameters:
    file (str): The path to the file being processed.
    threshold (float, optional): Threshold energy value to filter resonances. Defaults to None.
    """

    for res in resonances:
        if res.best_fit is not None:
            if threshold is None or res.threshold == threshold:
                th_dir = threshold_dir(project_dir, res.threshold)
                peak_fit(res.best_fit, f"{th_dir}[{res.index}]{res.energy:.8f}.png")
    if threshold is None:
        print(f"Plots saved to {project_dir}")
    else:
        print(f"Plots saved to {threshold_dir(project_dir, threshold)}")


def resonance_summary_grid(project_dir, resonances, resonance_index=None, open_files=None):
    """
    Create a grid plot for each resonance, showing all peaks associated with it.

    Parameters:
    project_dir (str): The project directory to save the summary plots.      # todo: colour codewords for darker colours (where subsequent ones are more different)
    resonances (list): List of Resonance objects containing associated peaks.  # todo: either allow emin,emax optional inputs, or determine energy range to only exclude few points but show peaks well.
    """

    norm = Normalize(vmin=1, vmax=4)  # Log scale for rel_SSR_per_point (10^-4 to 10^-1)
    cmap = plt.cm.get_cmap("RdYlGn")  # Gradient from red to green

    resonances_to_plot = (
        [resonances[resonance_index]] if resonance_index is not None and 0 <= resonance_index < len(resonances)
        else resonances
    )

    for res in resonances_to_plot:
        fig, axs = plt.subplots(nrows=((len(res.peaks)+1) // 3) + 1, ncols=3, figsize=(18, 6 * (((len(res.peaks)+1) // 3) + 1)))
        axs = axs.flatten()

        combined_ax = axs[0]
        combined_ax.set_title("Combined DOS Points for All Contributing Peaks")
        for idx, peak in enumerate(res.peaks):
            combined_ax.scatter(peak.energy_array, peak.dos_array, color=get_root_color(peak.root), s=3)
        combined_ax.set_xlabel("E (a.u.)")
        combined_ax.set_ylabel("DOS")
        combined_ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        for idx, peak in enumerate(sorted(res.peaks, key=lambda p: p.root)):
            ax = axs[idx+1]
            x_data = peak.energy_array
            y_data = peak.dos_array

            if not peak.is_descending:
                x_smooth = np.linspace(min(x_data), max(x_data), 1000)
                y_smooth = peak.get_smooth_lorentzian_curve(x_smooth)

                log_rel_ssr = -np.log10(max(peak.rel_ssr_per_point, 0) + 1)
                fit_color = cmap(norm(log_rel_ssr))

                ax.plot(x_smooth, y_smooth, color=fit_color)
                ax.set_title(f"Root {peak.root}, E = {peak.energy():.6f}, G = {peak.fit_Gamma:.6f}, Err = {peak.rel_ssr_per_point:.3e}")
            else:
                ax.set_title(f"Root {peak.root}, E = {peak.energy():.6f}, \n[!] Energy descending with growing gamma [!]", color='red')

            ax.scatter(x_data, y_data, edgecolor=get_root_color(peak.root), facecolor='white')
            ax.set_xlabel("Energy (a.u.)")
            ax.set_ylabel("DOS")
            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style='sci', axis='x', scilimits=(-2, 2))  # Use scientific notation if numbers are too large/small

            if peak == res.best_fit:
                annotation = "<Selected>"
                ax.text(0.5, 0.95, annotation, transform=ax.transAxes, fontsize=16, verticalalignment='top', horizontalalignment='center')

                # fit_text = (f"E = {peak.energy():.6f}\nGamma = {peak.fit_Gamma:.6f}\n"
            #             f"A = {peak.fit_A:.6f}\ny0 = {peak.fit_y0:.6f}\n"
            #             f"Rel SSR = {peak.rel_ssr_per_point:.3e}\n")
            # if len(peak.slope_energies):
            #     ax.scatter([peak.energy_array[i] for i in peak.slope_energies[:3]], [peak.dos_array[i] for i in peak.slope_energies[:3]], edgecolor="cyan")
            #     ax.scatter(peak.slope_energies[3:6], [0 for i in peak.slope_energies[3:6]], edgecolor="cyan")
            #     fit_text += f"{peak.slope_energies[-1]}"
            # if peak.warning:
            #     fit_text += f"\nWarning: {peak.warning}"
            # bbox_color = "green" if peak == res.best_fit else "white"
            # ax.text(0.30, 0.05, fit_text, transform=ax.transAxes, fontsize=10,
            #         verticalalignment='bottom', bbox=dict(boxstyle="round", facecolor="white", edgecolor=bbox_color, alpha=0.5))

        # Hide unused subplots
        for ax in axs[len(res.peaks)+1:]:
            ax.axis('off')

        plt.tight_layout()
        output_file = f"{threshold_dir(project_dir, res.threshold)}[{res.index}]{res.energy:.8f}.png"
        plt.savefig(output_file)
        plt.close()
        if resonance_index is not None:
            open_file(output_file, open_files)


def plot_all_resonance_peaks(data, resonances, output_file, emin=None, emax=None, clustering_output=None):
    """
    Plot all fitted peaks of all resonances in a single wide scatter plot (DOS vs. energy).

    Parameters:
    resonances (list): List of Resonance objects containing peaks.
    output_file (str): Path to save the generated plot.
    """
    plt.figure(figsize=(240, 8))  # Make the plot wide for better energy resolution

    # Plot the clustering amount if provided
    if clustering_output:
        energy_grid, clustering_array, peak_energies = clustering_output
        for peak_energy in peak_energies:
            plt.axvline(
                peak_energy,
                color="red",
                linestyle="--",
                linewidth=1,
                label=f"Detected Peak at {peak_energy:.4f}",
                alpha=0.5
            )
        plt.plot(
            energy_grid,
            (clustering_array)/20,
            color="blue",
            linewidth=1.5,
            label="Energy Clustering",
            alpha=0.6
        )

    for key in data.keys():
        if type(key) is str and key.startswith("rho_"):  # Identify DOS arrays
            root = int(key[4:])
            # indices = data[f"cleaned_{root}"][2:-2]
            rho = data[key]#[indices]
            energy = data[root][1:-1]#[indices]  # Corresponding energy array for the root
            color = get_root_color(root)
            plt.scatter(
                energy,
                np.log10(np.clip(rho, 0, None) + 1), #+1*int(key[4:]),
                # color="gray",
                color=color,
                alpha=0.5,
                s=5,
            )

    # for i, resonance in enumerate(resonances):
    #     # color = colors[i % num_colors]
    #     for peak in resonance.peaks:
    #         color = colors[peak.root % num_colors]
    #         plt.scatter(
    #             peak.energy_array,
    #             np.log10(peak.dos_array)+0.3*int(peak.root),
    #             color=color,
    #             label=f"Root {peak.root}",
    #             s=5,
    #             alpha=0.5
    #         )

    plt.xlabel("Energy (a.u.)")
    plt.ylabel("log10(DOS)")
    if emin is None:
        emin = -1
    if emax is None:
        emax = 0
    plt.xlim(emin, emax)
    plt.ylim(-0.2, 5)
    tick_positions = np.linspace(emin, emax, int(1/0.01))
    plt.xticks(tick_positions)
    plt.title("All Fitted Peaks: DOS vs Energy")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    open_file(output_file)



executor = ThreadPoolExecutor(max_workers=1)
clustering_future = None  # Will store the background future

def start_clustering_background_preparation(data, thresholds):
    global clustering_future
    clustering_future = executor.submit(prepare_plot_arrays, data, thresholds)

def get_plot_arrays():
    global clustering_future
    if clustering_future is not None:
        result = clustering_future.result()  # Waits here if still running
        return result

def prepare_plot_arrays(data, thresholds):
    """
    Prepares energy and rho arrays for clustering analysis, segmented by given thresholds.

    Parameters:
    - data (dict): Parsed data containing energy and DOS arrays.
    - thresholds (list): Energy thresholds defining sections for clustering.

    Returns:
    - dict: Dictionary containing full arrays and threshold-segmented arrays.
    """
    all_energies = []
    all_rhos = []
    line_data = {threshold: [] for threshold in thresholds}
    gamma = data["gamma"][1:-1]

    for key in data.keys():
        if isinstance(key, str) and key.startswith("rho_"):
            root = int(key[4:])
            rho = data[key]
            energy = data[root][1:-1]
            all_energies.append(energy)
            all_rhos.append(rho)

            prev_threshold = -np.inf
            for threshold in thresholds:
                mask = (energy >= prev_threshold) & (energy <= threshold)
                if np.any(mask):  # Only append non-empty segments
                    line_data[threshold].append((gamma[mask], energy[mask]))
                prev_threshold = threshold


    # Concatenate lists of arrays into single large NumPy arrays
    all_energies = np.concatenate(all_energies)
    all_rhos = np.concatenate(all_rhos)
    all_log10_rhos = np.log10(np.clip(all_rhos, 0, None) + 1)

    plot_arrays = {"all_energies": all_energies, "all_log10_rhos": all_log10_rhos, "line_data": line_data}

    prev_threshold = -np.inf
    for threshold in thresholds:
        mask = (all_energies >= prev_threshold) & (all_energies <= threshold)
        plot_arrays[f"energies_{threshold}"] = all_energies[mask]
        plot_arrays[f"log10_rhos_{threshold}"] = all_log10_rhos[mask]
        prev_threshold = threshold

    return plot_arrays

clustering_cache = {}

def resonance_partitions_with_clustering(data, resonances, emin, emax, output_file, open_files, threshold_above=0, manual_range = False):
    """
    Plot the partitioned sections of each root based on fitted peaks.

    Parameters:
    data (dict): Parsed data containing energy, gamma, and DOS arrays.
    fitted_peaks_by_root (dict): Dictionary of fitted peaks organized by root.
    output_file (str): The path where the plot will be saved.
    """


    plot_arrays = get_plot_arrays()  # prepared arrays trimmed for efficient plotting

    start = time.time()

    cache_key = (emin, emax)
    res_thr = [r for r in resonances if emin <= r.energy <= emax and r.best_fit is not None]
    res_thr.sort(key=lambda r: r.energy)

    if not manual_range and res_thr:
        best_fit = res_thr[0].best_fit
        emin = max(best_fit.energy() - 10 * best_fit.fit_Gamma, emin)

    end = time.time()
    print(f"Emin recalc: {end - start}")
    start = end

    fig = plt.figure(figsize=(21, 12))
    gs = GridSpec(1, 2, width_ratios=[16, 4], height_ratios=[9])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)

    end = time.time()
    print(f"Ax setup: {end - start}")
    start = end

    background_lines = plot_arrays["line_data"][threshold_above]
    for gamma_seg, energy_seg in background_lines:
        ax1.plot(gamma_seg, energy_seg, color="lightgray")

    end = time.time()
    print(f"Gray lines: {end - start}")
    start = end

    for res in resonances:
        show = res.should_be_shown()
        if emin < res.energy < emax:
            if show is not False:
                ax1.scatter(res.best_fit.gamma_array, res.best_fit.energy_array, color=get_root_color(res.index), s=5)
                for peak in res.peaks:
                    ax1.plot(peak.gamma_array, peak.energy_array, color=get_root_color(res.index), linewidth=3, alpha=0.2)
                    annotation_color = 'red' if peak.is_descending else 'black'  # "peaks" based on descending sections are marked red
                    vertical_offset = 0.0016 * (emax-emin)
                    if emin < peak.energy()+vertical_offset < emax:
                        ax1.text(peak.fit_gamma, peak.energy() + vertical_offset, f"{res.index}R{peak.root}", fontsize=8, ha='center', va='bottom', color=annotation_color, fontweight="bold" if peak==res.best_fit else "normal")

    end = time.time()
    print(f"resonance sections: {end - start}")
    start = end

    rotation = Affine2D().rotate_deg(90)  # Rotate rhs plot 90 degrees counterclockwise
    for res in resonances:
        if emin < res.energy < emax and res.best_fit is not None:
            show = res.should_be_shown()
            if show is False:
                continue
            annotation_color = 'black' if show is True else 'red'  # resonances marked (show = None) as based only on descending sections are marked red
            ax1.axhline(res.energy, color=get_root_color(res.index), linestyle="--", linewidth=1, alpha=0.2)
            ax2.axhline(res.energy, color=get_root_color(res.index), linestyle="--", linewidth=1, alpha=0.5)
            ax2.text(-4, res.energy, f"  [{res.index}] {res.energy:.6f}", ha="left", va="bottom", fontsize=8, color=annotation_color)


    end = time.time()
    print(f"Resonance lines: {end - start}")
    start = end

    # todo: same for panorama
    ax2.plot(plot_arrays[f"energies_{threshold_above}"], plot_arrays[f"log10_rhos_{threshold_above}"], '.', markersize=2, color="gray", transform=rotation + ax2.transData)  # todo: same for panorama

    end = time.time()
    print(f"Clustering scatter: {end - start}")
    start = end

    ax1.set_xlabel("gamma")
    ax2.set_xlabel("log(DOS)")
    ax1.set_ylabel("Energy (a.u.)")
    ax1.set_ylim(emin, emax)
    # ax1.set_title("Partitioned Sections of DOS by Resonance")
    ax2.set_xlim(-4, 0.1)
    ax2.tick_params(left=False, labelleft=False)
    # plt.ioff()  # Turn off interactive mode if accidentally enabled
    plt.subplots_adjust(wspace=0)
    plt.savefig(output_file, pil_kwargs={'compress_level': 1})
    plt.close()

    end = time.time()
    print(f"Savefig: {end - start}")
    start = end

    open_file(output_file, open_files)

    end = time.time()
    print(f"Open file: {end - start}")
    start = end


# Debug function
def plot_partitions(data, fitted_peaks_by_root, output_file, points):
    """
    Plot the partitioned sections of each root based on fitted peaks.

    Parameters:
    data (dict): Parsed data containing energy, gamma, and DOS arrays.
    fitted_peaks_by_root (dict): Dictionary of fitted peaks organized by root.
    output_file (str): The path where the plot will be saved.
    """
    plt.figure(figsize=(32, 24))

    for root in data.keys():
        if type(root) is int:
            plt.plot(data["gamma"], data[root], label=f"Root {root}", alpha=0.2)

    for root, peaks in fitted_peaks_by_root.items():
        for peak in peaks:
            plt.scatter(peak.gamma_array, peak.energy_array, label=f"Root {root} Section", s=5)

    plt.scatter([p[0] for p in points], [p[1] for p in points], s=7, facecolor="red", edgecolor="red")

            # vertical_offset = -0.0015 if peak.is_left_half else 0.0015
            # plt.text(
            #     peak.fit_gamma,  # X-coordinate
            #     peak.energy() + vertical_offset,  # Y-coordinate with offset
            #     # f"E={peak.energy():.3f}\nρ={peak.fit_Gamma:.3f}\n{peak.fit_A:.3f}, {peak.fit_y0:.3f}",  # Annotation text
            #     # f"{peak.fit_Gamma:.1e}, {peak.fit_y0:.3f}, {peak.max_dos:.3f}",  # Annotation text
            #     f"{peak.root}",  # Annotation text
            #     fontsize=8,
            #     ha='center',  # Horizontal alignment
            #     va='bottom' if peak.is_left_half else 'top',  # Vertical alignment
            #     color='black'
            # )

    plt.xlabel("Gamma")
    plt.ylabel("Energy (a.u.)")
    plt.ylim(-0.7, -0.5)
    plt.title("Partitioned Sections of DOS by Root")
    plt.savefig(output_file)
    plt.close()
    open_file(output_file)


