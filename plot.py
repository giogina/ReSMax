import os
import subprocess
import platform
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter


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
    open_plot(file)


def peak_fit(dos_peak, file):
    """
    Plot the fitted peak for a DOSpeak object and save the plot.

    Parameters:
    dos_peak (DOSpeak): The DOSpeak object containing the fitted data.
    file (str): The path where the plot should be saved.
    """
    x_data = dos_peak.energy_array
    y_data = dos_peak.dos_array
    x_smooth = np.linspace(min(min(x_data), dos_peak.fit_E - dos_peak.fit_Gamma/2), max(x_data), 1000) if dos_peak.discontinuity_treatment == "fit" else np.linspace(min(x_data), max(x_data), 1000)
    y_smooth = dos_peak.get_smooth_lorentzian_curve(x_smooth)
    plt.figure(figsize=(12, 8))
    plt.scatter(x_data, y_data, edgecolor='black', facecolor='none')
    plt.plot(x_smooth, y_smooth, 'r-')
    plt.xlabel("Energy (a.u.)")
    plt.ylabel("DOS")
    plt.savefig(file)
    plt.close()
    with open(f"{file[:-5]}.txt",'w') as plot_data_file:
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
    plt.figure(figsize=(16, 12))
    min_E = 0 if from_e is None else from_e
    max_E = -100 if to_e is None else to_e

    for key, values in data.items():
        if key != "gamma":
            plt.plot(data["gamma"], values, label=key)
            if from_e is None:
                min_E = min(min_E, min(values))
            if to_e is None:
                max_E = max(max_E, max(values))

    plt.xlabel("γ")
    plt.ylabel("Energy (a.u.)")
    plt.ylim(bottom=min_E - 0.03, top=max_E + 0.03)
    plt.locator_params(axis='x', nbins=20)
    plt.locator_params(axis='y', nbins=20)
    plt.minorticks_on()

    plt.savefig(plot_file)
    plt.close()
    open_plot(plot_file)
    return min_E, max_E


def open_plot(file):
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


def resonance_fits(project_dir, resonances, threshold=None):
    """
    Plot the resonance fits and save them in the appropriate directories.

    Parameters:
    file (str): The path to the file being processed.
    threshold (float, optional): Threshold energy value to filter resonances. Defaults to None.
    """

    sep = '\\' if '\\' in project_dir else '/'
    plots_dir = f"{project_dir}resonance_plots{sep}"
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    if threshold is None:
        for t in set([res.threshold for res in resonances]):
            if not os.path.exists(f"{plots_dir}{t}{sep}"):
                os.mkdir(f"{plots_dir}{t}{sep}")
    else:
        if not os.path.exists(f"{plots_dir}{threshold}{sep}"):
            os.mkdir(f"{plots_dir}{threshold}{sep}")
    for res in resonances:
        if threshold is None or res.threshold == threshold:
            peak_fit(res.best_fit, f"{plots_dir}{res.threshold}{sep}{res.energy}.png")
    if threshold is None:
        print(f"Plots saved to {plots_dir}")
    else:
        print(f"Plots saved to {plots_dir}{threshold}{sep}")


def resonance_summary_grid(project_dir, resonances):
    """
    Create a grid plot for each resonance, showing all peaks associated with it.

    Parameters:
    project_dir (str): The project directory to save the summary plots.
    resonances (list): List of Resonance objects containing associated peaks.
    """
    sep = '\\' if '\\' in project_dir else '/'
    summary_dir = f"{project_dir}resonance_plot_grids{sep}"
    if not os.path.exists(summary_dir):
        os.mkdir(summary_dir)

    norm = Normalize(vmin=1, vmax=4)  # Log scale for rel_SSR_per_point (10^-4 to 10^-1)
    cmap = plt.cm.get_cmap("RdYlGn")  # Gradient from red to green

    for res in resonances:
        fig, axs = plt.subplots(nrows=(len(res.peaks) // 3) + 1, ncols=3, figsize=(18, 6 * ((len(res.peaks) // 3) + 1)))
        axs = axs.flatten()
        for idx, peak in enumerate(res.peaks):
            ax = axs[idx]
            x_data = peak.energy_array
            y_data = peak.dos_array
            x_smooth = np.linspace(min(x_data), max(x_data), 1000)
            y_smooth = peak.get_smooth_lorentzian_curve(x_smooth)

            # Determine fit curve color based on log(rel_SSR_per_point)
            log_rel_ssr = -np.log10(peak.rel_ssr_per_point)
            fit_color = cmap(norm(log_rel_ssr))

            # Data points color: red if warning exists, otherwise black
            points_color = 'red' if peak.warning else 'black'

            # Plotting
            ax.scatter(x_data, y_data, edgecolor=points_color, facecolor='white')
            ax.plot(x_smooth, y_smooth, color=fit_color)
            ax.set_title(f"Root {peak.root}, E_fit = {peak.fit_E:.5f}")
            ax.set_xlabel("Energy (a.u.)")
            ax.set_ylabel("DOS")
            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style='sci', axis='x',
                                scilimits=(-2, 2))  # Use scientific notation if numbers are too large/small

            # Annotate fit parameters and warnings
            fit_text = (f"E = {peak.fit_E:.6f}\nGamma = {peak.fit_Gamma:.6f}\n"
                        f"A = {peak.fit_A:.6f}\ny0 = {peak.fit_y0:.6f}\n"
                        f"Rel SSR = {peak.rel_ssr_per_point:.3e}")
            # if peak.warning:
            #     fit_text += f"\nWarning: {peak.warning}"
            bbox_color = "green" if peak == res.best_fit else "white"
            ax.text(0.30, 0.05, fit_text, transform=ax.transAxes, fontsize=10,
                    verticalalignment='bottom', bbox=dict(boxstyle="round", facecolor="white", edgecolor=bbox_color, alpha=0.5))

        # Hide unused subplots
        for ax in axs[len(res.peaks):]:
            ax.axis('off')

        plt.tight_layout()
        plt.savefig(f"{summary_dir}{res.energy:.5f}.png")
        plt.close()
