import os
import subprocess
import platform
import matplotlib.pyplot as plt
import numpy as np


def plot_DOS(data, root, file, fitted_peaks_by_root=None):
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


def overview(data, plot_file, from_e=None, to_e=0):
    plt.figure(figsize=(16, 12))
    min_E = 0 if from_e is None else from_e
    max_E = -100

    for key, values in data.items():
        if key != "gamma":
            plt.plot(data["gamma"], values, label=key)
            if from_e is None:
                min_E = min(min_E, min(values))
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
