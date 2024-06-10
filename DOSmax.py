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
                            max_e = -z * z / 2. / 36  # should never happen - limit to 6 threshold values
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a file.")
    parser.add_argument('-f', '--file', type=str, required=True, help="Path to the file")
    args = parser.parse_args()
    main(args.file)
