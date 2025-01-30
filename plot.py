import os
import subprocess
import platform
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import Affine2D


def get_root_color(rootnr: int):
    colors = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan",
              "navy", "teal", "gold", "lime", "maroon", "yellow", "magenta", "turquoise", "indigo"]
    return colors[rootnr % len(colors)]


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
    xmin, xmax = dos_peak.fit_E-8*dos_peak.fit_Gamma, dos_peak.fit_E+8*dos_peak.fit_Gamma
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
    with open_file(f"{file[:-5]}.txt", 'w') as plot_data_file:
        plot_data_file.write(np.array2string(x_data, separator=" ", max_line_width=np.inf) + "\r\n")
        plot_data_file.write(np.array2string(y_data, separator=" ", max_line_width=np.inf) + "\r\n")
        plot_data_file.write(np.array2string(x_smooth, separator=" ", max_line_width=np.inf) + "\r\n")
        plot_data_file.write(np.array2string(y_smooth, separator=" ", max_line_width=np.inf) + "\r\n")


def overview(data, plot_file, from_e=None, to_e=None):
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
            plt.plot(data["gamma"], values, label=key)  # plot / scatter
            # if len(data[f"cleaned_{key}"]):
            #     plt.scatter(data["gamma"][data[f"cleaned_{key}"]], values[data[f"cleaned_{key}"]], label=key, s=5, color=get_root_color(key))
            # plt.scatter(data["gamma"], values, label=key, s=5, color=get_root_color(key))
            if from_e is None:
                min_E = min(min_E, min(values))
            if to_e is None:
                max_E = max(max_E, max(values))

    plt.xlabel("γ")
    plt.ylabel("Energy (a.u.)")
    plt.ylim(bottom=min_E - 0.03, top=max_E + 0.03)
    # plt.ylim(bottom=min_E - 0.0, top=max_E + 0.0)
    # plt.ylim(bottom=-0.38, top=-0.22)
    plt.locator_params(axis='x', nbins=20)
    plt.locator_params(axis='y', nbins=20)
    plt.minorticks_on()

    plt.savefig(plot_file)
    plt.close()
    open_file(plot_file)
    return min_E, max_E


def open_file(file):
    """
    Open a file using the default application based on the operating system.

    Parameters:
    file (str): The file path to open.
    """
    if platform.system() == 'Windows':
        os.startfile(file)
    elif platform.system() == 'Darwin':  # macOS
        subprocess.call(('open', file))
    else:  # Linux and other Unix-like systems
        subprocess.call(('xdg-open', file))


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


def resonance_summary_grid(project_dir, resonances, resonance_index=None):
    """
    Create a grid plot for each resonance, showing all peaks associated with it.

    Parameters:
    project_dir (str): The project directory to save the summary plots.
    resonances (list): List of Resonance objects containing associated peaks.
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
            x_smooth = np.linspace(min(x_data), max(x_data), 1000)
            y_smooth = peak.get_smooth_lorentzian_curve(x_smooth)

            log_rel_ssr = -np.log10(peak.rel_ssr_per_point)
            fit_color = cmap(norm(log_rel_ssr))

            ax.scatter(x_data, y_data, edgecolor=get_root_color(peak.root), facecolor='white')
            ax.plot(x_smooth, y_smooth, color=fit_color)
            ax.set_title(f"Root {peak.root}, E = {peak.fit_E:.6f}, G = {peak.fit_Gamma:.6f}, Err = {peak.rel_ssr_per_point:.3e}")
            ax.set_xlabel("Energy (a.u.)")
            ax.set_ylabel("DOS")
            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style='sci', axis='x', scilimits=(-2, 2))  # Use scientific notation if numbers are too large/small

            # fit_text = (f"E = {peak.fit_E:.6f}\nGamma = {peak.fit_Gamma:.6f}\n"
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
            open_file(output_file)


def plot_all_resonance_peaks(data, resonances, output_file, clustering_output=None):
    """
    Plot all fitted peaks of all resonances in a single wide scatter plot (DOS vs. energy).

    Parameters:
    resonances (list): List of Resonance objects containing peaks.
    output_file (str): Path to save the generated plot.
    """
    plt.figure(figsize=(240, 8))  # Make the plot wide for better energy resolution

    colors = ["blue", "purple", "red", "yellow", "green", "cyan"]
    num_colors = len(colors)

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
            # indices = data[f"cleaned_{root}"][2:-2]  # todo: peaks are skewed now; this seems wrong?
            rho = data[key]#[indices]
            energy = data[root][1:-1]#[indices]  # Corresponding energy array for the root
            color = get_root_color(root)
            plt.scatter(
                energy,
                np.log10(rho+1), #+1*int(key[4:]),
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
    # emin, emax = -0.53, -0.50
    emin, emax = -1, 0
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


def plot_resonance_partitions_with_clustering(data, resonances, emin, emax, output_file):
    """
    Plot the partitioned sections of each root based on fitted peaks.

    Parameters:
    data (dict): Parsed data containing energy, gamma, and DOS arrays.
    fitted_peaks_by_root (dict): Dictionary of fitted peaks organized by root.
    output_file (str): The path where the plot will be saved.
    """

    res_thr = [r for r in resonances if emin <= r.energy <= emax and r.best_fit is not None]
    res_thr.sort(key=lambda r: r.energy)
    emin = max(res_thr[0].best_fit.fit_E - 10 * res_thr[0].best_fit.fit_Gamma, emin)

    fig = plt.figure(figsize=(21, 12))
    gs = GridSpec(1, 2, width_ratios=[16, 4], height_ratios=[9])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)

    for root in data.keys():
        if type(root) is int:
            ax1.plot(data["gamma"], data[root], color="gray", alpha=0.2)

    for res in resonances:
        if emin < res.energy < emax:
            if res.best_fit is not None:
                ax1.scatter(res.best_fit.gamma_array, res.best_fit.energy_array, color=get_root_color(res.index), s=5)
            for peak in res.peaks:
                ax1.scatter(peak.gamma_array, peak.energy_array, color=get_root_color(res.index), s=5, alpha=0.1)
                vertical_offset = 0.0005
                if emin < peak.fit_E+vertical_offset < emax:
                    ax1.text(peak.fit_gamma, peak.fit_E + vertical_offset, f"{res.index}R{peak.root}", fontsize=8, ha='center', va='bottom', color='black')

    rotation = Affine2D().rotate_deg(90)  # Rotate rhs plot 90 degrees counterclockwise
    for res in resonances:
        if emin < res.energy < emax and res.best_fit is not None:
            ax1.axhline(res.energy, color=get_root_color(res.index), linestyle="--", linewidth=1, alpha=0.2)
            ax2.axhline(res.energy, color=get_root_color(res.index), linestyle="--", linewidth=1, alpha=0.5)
            ax2.text(-4, res.energy, f"  [{res.index}] {res.energy:.6f}", ha="left", va="bottom", fontsize=8)

    for key in data.keys():
        if type(key) is str and key.startswith("rho_"):  # Identify DOS arrays
            root = int(key[4:])
            rho = data[key]
            energy = data[root][1:-1]
            ax2.scatter(energy, np.log10(rho + 1), color="gray", alpha=0.5, s=5, transform=rotation + ax2.transData)

    ax1.set_xlabel("gamma")
    ax2.set_xlabel("log(DOS)")
    ax1.set_ylabel("Energy (a.u.)")
    ax1.set_ylim(emin, emax)
    # ax1.set_title("Partitioned Sections of DOS by Resonance")
    ax2.set_xlim(-4, 0.1)
    ax2.tick_params(left=False, labelleft=False)
    plt.subplots_adjust(wspace=0)
    plt.savefig(output_file)
    plt.close()
    open_file(output_file)

