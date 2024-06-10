import os
from scipy import signal
import argparse
import numpy as np



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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a file.")
    parser.add_argument('-f', '--file', type=str, required=True, help="Path to the file")
    args = parser.parse_args()
    main(args.file)
