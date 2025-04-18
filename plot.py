import gc
import os
import subprocess
import platform

import numpy as np
import psutil
from concurrent.futures import ThreadPoolExecutor

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, hsv_to_rgb
from matplotlib.gridspec import GridSpec
from matplotlib.text import Text
from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MultipleLocator, FormatStrFormatter, MaxNLocator, \
    AutoLocator, LogLocator, NullFormatter, FixedLocator
from matplotlib.transforms import Affine2D


def get_root_color(rootnr: int, alpha=1.):
    saturation = 1
    value = 0.7
    hue = (rootnr * 222.5+15) % 360
    hue /= 360.0
    rgb = hsv_to_rgb((hue, saturation, value))
    if alpha != 1:
        rgb = (1 - alpha) * np.array([1, 1, 1]) + alpha * rgb
    return rgb

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
    plt.close('all')
    open_file(file)


def peak_fit(dos_peak, file):
    """
    Plot the fitted peak for a DOSpeak object and save the plot.

    Parameters:
    dos_peak (DOSpeak): The DOSpeak object containing the fitted data.
    file (str): The path where the plot should be saved.
    """

    plt.rcParams.update({
        "font.size": 14,  # default font size
        "axes.titlesize": 16,  # title font
        "axes.labelsize": 14,  # x/y label font
        "xtick.labelsize": 14,  # x-tick font
        "ytick.labelsize": 14,  # y-tick font
        "legend.fontsize": 14,
    })

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
    plt.minorticks_on()
    # plt.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.xlim(xmin, xmax)
    plt.savefig(file)
    plt.close('all')
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
    # plt.figure(figsize=(32, 24))
    plt.figure(figsize=(16, 12))
    min_E = 0 if from_e is None else from_e
    max_E = -100 if to_e is None else to_e

    for key, values in data.items():
        if type(key) is int:
            plt.plot(data["gamma"], values, label=key, color=get_root_color(key))  # plot / scatter
            if from_e is None:
                min_E = min(min_E, min(values))
            if to_e is None:
                max_E = max(max_E, max(values))

    labelSize = 20
    plt.xlabel("γ", fontsize=labelSize, labelpad=int(1.0*labelSize))
    plt.ylabel("Energy (a.u.)", fontsize=labelSize, labelpad=int(1.0*labelSize))
    plt.tick_params(axis='both', which='major', labelsize=labelSize)
    plt.minorticks_on()
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
    plt.ylim(bottom=min_E - margin, top=max_E + margin)
    plt.tick_params(axis='both', which='major', length=10, width=2, pad=10)  # Major ticks larger
    plt.tick_params(axis='both', which='minor', length=5, width=1, pad=10)  # Minor ticks slightly smaller
    plt.locator_params(axis='x', nbins=20)
    plt.locator_params(axis='y', nbins=20)
    ax = plt.gca()  # Get current axis
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.tight_layout()

    plt.savefig(plot_file)
    plt.clf()  # Clear figure
    plt.cla()
    plt.close('all')  # Extra call just in case
    gc.collect()
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
        if res.best_fit is not None and not res.best_fit.is_descending:
            if threshold is None or res.threshold == threshold:
                th_dir = threshold_dir(project_dir, res.threshold)
                peak_fit(res.best_fit, f"{th_dir}[{res.index+1}]{res.energy:.8f}.png")
    if threshold is None:
        print(f"Plots saved to {project_dir}")
    else:
        print(f"Plots saved to {threshold_dir(project_dir, threshold)}")


def resonance_summary_grid(project_dir, resonances, resonance_index=None, open_files=None):
    """
    Create a grid plot for each resonance, showing all peaks associated with it.

    Parameters:
    project_dir (str): The project directory to save the summary plots
    resonances (list): List of Resonance objects containing associated peaks.
    """

    # norm = Normalize(vmin=1, vmax=4)  # Log scale for rel_SSR_per_point (10^-4 to 10^-1)
    # cmap = plt.cm.get_cmap("RdYlGn")  # Gradient from red to green

    resonances_to_plot = (
        [resonances[resonance_index]] if resonance_index is not None and 0 <= resonance_index < len(resonances)
        else resonances
    )

    grid_width = 2

    plt.rcParams.update({
        "font.size": 16,  # default font size
        "axes.titlesize": 18,  # title font
        "axes.labelsize": 16,  # x/y label font
        "xtick.labelsize": 16,  # x-tick font
        "ytick.labelsize": 16,  # y-tick font
        "legend.fontsize": 16,
    })

    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, 3))
    formatter.set_useOffset(False)  # disables +offset shifting

    for res in resonances_to_plot:
        nrrows = ((len(res.peaks)+1) // grid_width) + 1
        fig, axs = plt.subplots(nrows=nrrows, ncols=grid_width, figsize=(grid_width*6, 6*nrrows))
        axs = axs.flatten()

        combined_ax = axs[0]
        combined_ax.set_title("Combined DOS Points\nfor All Contributing Peaks")
        threshold = 0.005 * max([np.max(p.dos_array, 0) for p in res.peaks])
        for idx, peak in enumerate(res.peaks):
            valid_indices = np.where(np.array(peak.dos_array) >= threshold)[0]
            combined_ax.plot(peak.energy_array[valid_indices], peak.dos_array[valid_indices], '.', markersize=2, color=get_root_color(peak.root))
        combined_ax.set_xlabel("Energy (a.u.)")
        combined_ax.set_ylabel("DOS")
        # combined_ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        combined_ax.xaxis.set_major_formatter(formatter)
        combined_ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))
        combined_ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        combined_ax.ticklabel_format(style='sci', axis='y', scilimits=(-1, 1))  # Use scientific notation if numbers are too large/small
        combined_ax.yaxis.get_offset_text().set_horizontalalignment('right')

        majorticks = combined_ax.get_xticks()
        step = np.diff(majorticks).min()
        minor_step = 10 ** (np.ceil(np.log10(step)) - 1)
        xmin, xmax = combined_ax.get_xlim()
        minor_ticks = np.arange(np.floor(xmin / minor_step) * minor_step, np.ceil(xmax / minor_step) * minor_step, minor_step)
        # print(step, minor_step, xmin, xmax, np.floor(xmin / minor_step) * minor_step, np.ceil(xmax / minor_step) * minor_step)

        combined_ax.xaxis.set_minor_locator(FixedLocator(minor_ticks))
        combined_ax.xaxis.set_minor_formatter(NullFormatter())
        combined_ax.tick_params(axis='both', which='major', length=6, width=1.5)
        combined_ax.tick_params(axis='both', which='minor', length=4, width=1.0)

        sorted_peaks = sorted(res.peaks, key=lambda p: p.root)
        # log_rel_ssrs = np.array([-np.log10(max(peak.rel_ssr_per_point, 0) + 1) for peak in sorted_peaks])
        # fit_colors = cmap(norm(log_rel_ssrs))
        fit_color = "black"
        for idx, peak in enumerate(sorted_peaks):
            ax = axs[idx+1]
            x_data = peak.energy_array
            y_data = peak.dos_array

            if not peak.is_descending:
                x_smooth = np.linspace(min(x_data), max(x_data), 1000)
                y_smooth = peak.get_smooth_lorentzian_curve(x_smooth)
                ax.plot(x_smooth, y_smooth, color=fit_color)
                # y_smooth_guess = peak.get_smooth_guess_lorentzian_curve(x_smooth)  # debug plots for checking initial guesses - looking good!
                # ax.plot(x_smooth, y_smooth_guess, color="blue")
                ax.set_title(f"Root {peak.root}, E = {peak.energy():.6f}\nG = {peak.fit_Gamma:.6f}, Err = {peak.rel_ssr_per_point:.3e}")
            else:
                ax.set_title(f"Root {peak.root}, E = {peak.energy():.6f}\n[!] Energy descending [!]\n with growing gamma", color='red')

            ax.scatter(x_data, y_data, edgecolor=get_root_color(peak.root), facecolor='white')
            ax.set_xlabel("Energy (a.u.)")
            ax.set_ylabel("DOS")

            ax.xaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))

            xmin, xmax = ax.get_xlim()
            minor_ticks = np.arange(np.floor(xmin / minor_step) * minor_step, np.ceil(xmax / minor_step) * minor_step, minor_step)
            ax.xaxis.set_minor_locator(FixedLocator(minor_ticks))
            ax.xaxis.set_minor_formatter(NullFormatter())

            ax.tick_params(axis='both', which='major', length=6, width=1.5)
            ax.tick_params(axis='both', which='minor', length=4, width=1.0)

            # ax.yaxis.set_major_formatter(FormatStrFormatter('%.0e'))
            # ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            # ax.ticklabel_format(style='sci', axis='x', scilimits=(-1, 1))  # Use scientific notation if numbers are too large/small
            ax.ticklabel_format(style='sci', axis='y', scilimits=(-1, 1))  # Use scientific notation if numbers are too large/small
            ax.yaxis.get_offset_text().set_horizontalalignment('right')

            if peak == res.best_fit:
                annotation = "<Selected>"
                ax.text(0.1, 0.95, annotation, transform=ax.transAxes, fontsize=18, verticalalignment='top', horizontalalignment='left')

        # Hide unused subplots
        for ax in axs[len(res.peaks)+1:]:
            ax.axis('off')

        plt.tight_layout()
        output_file = f"{threshold_dir(project_dir, res.threshold)}[{res.index+1}]{res.energy:.8f}.png"
        fig.savefig(output_file)
        plt.close('all')
        if resonance_index is not None:
            open_file(output_file, open_files)


def plot_all_resonance_peaks(data, resonances, output_file, emin=None, emax=None, clustering_output=None):
    """
    Plot all DOS values in a single wide scatter plot (DOS vs. energy).

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
            plt.plot(
                energy,
                np.log10(np.clip(rho, 0, None) + 1), #+1*int(key[4:]),
                '.',
                markersize=2,
                color=color,
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
    plt.close('all')
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


def resonance_partitions_with_clustering(data, resonances, emin, emax, output_file, open_files, threshold_above=0, manual_range = False):
    """
    Plot the partitioned sections of each root based on fitted peaks.

    Parameters:
    data (dict): Parsed data containing energy, gamma, and DOS arrays.
    fitted_peaks_by_root (dict): Dictionary of fitted peaks organized by root.
    output_file (str): The path where the plot will be saved.
    """


    plot_arrays = get_plot_arrays()  # prepared arrays trimmed for efficient plotting
    executor.shutdown(wait=True)
    # print("plot arrays: ", plot_arrays)
    res_thr = [r for r in resonances if emin <= r.energy <= emax]
    res_thr.sort(key=lambda r: r.energy)

    if not manual_range and res_thr:
        best_fit = res_thr[0].best_fit
        if best_fit is not None and best_fit.fit_Gamma is not None:
            emin = max(best_fit.energy() - 10 * best_fit.fit_Gamma, emin)

    # fig = plt.figure(figsize=(21, 12))
    fig = plt.figure(figsize=(15, 5))
    gs = GridSpec(1, 2, width_ratios=[16, 4], height_ratios=[9])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)

    background_lines = plot_arrays["line_data"][threshold_above]
    for gamma_seg, energy_seg in background_lines:
        ax1.plot(gamma_seg, energy_seg, color="lightgray")

    texts_to_add = []
    for res in res_thr:
        show = res.should_be_shown()
        res_gammas = []
        res_es = []
        for peak in res.peaks:
            if len(peak.gamma_array):
                res_gammas.append(peak.gamma_array)
                res_gammas.append([np.nan])
                res_es.append(peak.energy_array)
                res_es.append([np.nan])
                annotation_color = 'red' if peak.is_descending else 'black'  # "peaks" based on descending sections are marked red
                vertical_offset = 0.0016 * (emax-emin)
                if emin < peak.energy()+vertical_offset < emax:
                    # ax1.text(peak.fit_gamma, peak.energy() + vertical_offset, f"{res.index}R{peak.root}", fontsize=8, ha='center', va='bottom', color=annotation_color, fontweight="bold" if peak==res.best_fit else "normal")
                    texts_to_add.append(Text(x=peak.fit_gamma, y=peak.energy() + vertical_offset, text=f"{res.index+1}R{peak.root}", fontsize=8, ha='center', va='bottom', color=annotation_color, fontweight="bold" if peak==res.best_fit and res.should_be_shown() else "normal"))

        res_gamma_all = np.concatenate(res_gammas) if res_gammas else np.array([])
        res_energy_all = np.concatenate(res_es) if res_es else np.array([])
        ax1.plot(res_gamma_all, res_energy_all, color=get_root_color(res.index, alpha=0.2), linewidth=3, solid_capstyle="round")

        if show is not False:
            ax1.plot(res.best_fit.gamma_array, res.best_fit.energy_array, color=get_root_color(res.index), linewidth=3, solid_capstyle="round")

    for text_obj in texts_to_add:
        ax1.add_artist(text_obj)

    rotation = Affine2D().rotate_deg(90)  # Rotate rhs plot 90 degrees counterclockwise
    ax2.plot(plot_arrays[f"energies_{threshold_above}"], plot_arrays[f"log10_rhos_{threshold_above}"], '.', markersize=2, color="gray", transform=rotation + ax2.transData)

    text_height = 0.03*(emax-emin)
    prev_res = None
    shown_resonances = [r for r in resonances if emin < r.energy < emax and r.best_fit is not None and (r.should_be_shown() is not False)]
    for i, res in enumerate(shown_resonances):
        show = res.should_be_shown()
        annotation_color = 'black' if show is True else 'red'  # resonances marked (show = None) as based only on descending sections are marked red
        ax1.axhline(res.energy, color=get_root_color(res.index), linestyle="--", linewidth=1, alpha=0.2)
        ax2.axhline(res.energy, color=get_root_color(res.index, alpha=0.5), linestyle="--", linewidth=1)

        next_res = shown_resonances[i+1] if (i+1<len(shown_resonances)) else None
        if (next_res is not None) and (next_res.energy-res.energy < text_height) and ((prev_res is None) or (res.energy - prev_res.energy > 2*text_height)):
            ax2.text(-4, res.energy-0.0*text_height, f"  [{res.index+1}] {res.energy:.6f}", ha="left", va="top", fontsize=8, color=annotation_color)
        else:
            ax2.text(-4, res.energy+0.0*text_height, f"  [{res.index+1}] {res.energy:.6f}", ha="left", va="bottom", fontsize=8, color=annotation_color)
        prev_res = res


    for artist in ax1.get_children():
        if hasattr(artist, 'set_rasterized'):
            artist.set_rasterized(True)
    for artist in ax2.get_children():
        if hasattr(artist, 'set_rasterized'):
            artist.set_rasterized(True)

    ax1.set_xlabel("γ")
    ax2.set_xlabel("log(DOS)")
    ax1.set_ylabel("Energy (a.u.)")
    ax1.set_ylim(emin, emax)

    # plt.tick_params(axis='both', which='major', labelsize=12, length=12, width=2, pad=15)  # Major tick adjustments
    # plt.tick_params(axis='both', which='minor', length=6, width=1.5)
    # ax1.set_title("Partitioned Sections of DOS by Resonance")
    ax2.set_xlim(-4, 0.1)
    ax2.tick_params(left=False, labelleft=False)
    plt.minorticks_on()
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    plt.subplots_adjust(wspace=0)
    plt.subplots_adjust(bottom=0.15)

    fig.savefig(output_file, pil_kwargs={'compress_level': 1})
    plt.close('all')
    open_file(output_file, open_files)


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
    # plt.ylim(-0.52, -0.515)
    plt.title("Partitioned Sections of DOS by Root")
    plt.savefig(output_file)
    plt.close('all')
    open_file(output_file)


