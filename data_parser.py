import re
import os
import numpy as np


float_pattern = re.compile(r'^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$')


def parse(file: str):
    """
    Parse the input file based on its extension and return the data structure.

    Parameters:
    file (str): The path to the file.

    Returns:
    dict: Parsed data from the file.
    """
    if file.endswith(".ou"):
        return parse_ou(file)
    elif file.endswith(".dat"):
        return parse_dat(file)
    elif file.endswith(".dal"):
        return parse_dal(file)
    else:
        print("Unrecognized data type. Accepted: .ou (gamma array followed by E(gamma, root) array for each root) or .dat (arrays of gamma, E(gamma, root) for each gamma)")
        return None


def parse_dat(file):
    """
    Parse a .dat file and return the data in a dictionary.
    Data shape: Tabular, columns = [gamma, E[root 1](gamma), E[root 2](gamma), ...]

    Parameters:
    file (str): The path to the .dat file.

    Returns:
    dict: Parsed data from the .dat file.
    """
    res = {"gamma": []}
    nr_data_points = -1
    with open(file) as f:
        lines = f.readlines()
    if lines is not None:
        for root, e in enumerate(lines[0].split()[1:]):
            res[root + 1] = []  # initiate root energy arrays
        for line in lines:
            vals = line.split()
            res["gamma"].append(float(vals[0]))
            for root, e in enumerate(vals[1:]):
                res[root + 1].append(float(e))
    nr_data_points = len(res["gamma"])
    for root in res.keys():
        if len(res[root]) == nr_data_points:
            res[root] = np.array(res[root])
        else:
            print(f"Wrong number of data points for root #{root}: {len(res[root])}, should be {nr_data_points}.")
            res.pop(root)
    return res


def parse_dal(file):
    """
    Parse a .dal file and return the data in a dictionary.
    Data shape: One number per line, first gamma, then E[root 1](gamma), E[root 2](gamma), .... Double new line, then next block.

    Parameters:
    file (str): The path to the .dal file.

    Returns:
    dict: Parsed data from the .dal file.
    """
    res = {"gamma": []}
    with open(file) as f:
        lines = f.readlines()
    if lines is not None:
        vals = []
        for line in lines:
            if len(line.strip()) == 0:
                if len(vals):
                    gamma = float(vals[0])
                    res["gamma"].append(gamma)
                    for root, e in enumerate(vals[1:]):
                        if not root+1 in res.keys():
                            res[root + 1] = {}  # initiate root energy dict
                        res[root + 1][gamma] = float(e)
                vals = []
            else:
                vals.append(float(line.strip()))

    res["gamma"] = np.array(sorted(res["gamma"]))
    nr_data_points = len(res["gamma"])
    roots = res.keys()
    for root in roots:
        if type(root) is int:
            root_vals = [res[root][g] for g in sorted(res[root].keys())]
            if len(root_vals) == nr_data_points:
                res[root] = np.array(root_vals)
            else:
                print(f"Wrong number of data points for root #{root}: {len(root_vals)}, should be {nr_data_points}.")
                res.pop(root)
    return res


def parse_ou(file):
    """
    Parse a .ou file and return the data in a dictionary.

    Parameters:
    file (str): The path to the .ou file.

    Returns:
    dict: Parsed data from the .ou file.
    """
    res = dict()
    nr_data_points = -1
    with open(file) as f:
        lines = f.readlines()
    if lines is not None:
        data_type = "gamma"
        data = list()
        for line in lines:
            if is_float(line):
                data.append(float(line.strip()))
            else:
                if len(data):
                    if data_type == "gamma":
                        nr_data_points = len(data)
                        res[data_type] = np.array(data.copy())
                        data_type = 1
                    else:
                        if len(data) == nr_data_points:
                            res[data_type] = np.array(data.copy())
                        data_type += 1
                    data = list()
    return res


def is_float(value: str):
    """
    Check if a string represents a floating-point number.

    Parameters:
    value (str): The string to check.

    Returns:
    bool: True if the string represents a float, otherwise False.
    """
    return bool(float_pattern.match(value.strip()))


def project_directory(file: str):
    """
    Create or retrieve the project directory based on the file path.

    Parameters:
    file (str): The file path.

    Returns:
    str: The project directory path.
    """
    sep = '\\' if '\\' in file else '/'
    project_dir = os.path.splitext(file)[0] + sep
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)
    return project_dir
