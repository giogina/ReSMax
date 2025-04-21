# DOSmax
**Efficient, automated resonance detection using the Stabilization Method.**

**DOSmax** is a Python tool for identifying resonance states in quantum few-body systems, such as helium-like atoms, using the Stabilization Method. This technique analyzes how energy eigenvalues evolve with a systematically varied basis set parameter, revealing metastable bound states as well as resonance states, which appear as characteristic plateaus in the eigenvalue spectrum.

**DOSmax** automates the full workflow: it computes the density of states (DOS) for each root, segments and fits DOS peaks using Lorentzian profiles, groups them into resonances, and extracts resonance parameters (energy and width) from the best-fitting DOS peaks — all with minimal manual input. Results can be interactively inspected and refined via plots and command-line controls, allowing for careful manual validation.

The tool is designed for reproducibility and efficient large-scale analysis of stabilization data, helping researchers extract physical resonance parameters quickly and consistently.

## Features
- Automated resonance detection from stabilization diagrams with minimal user input.
- Supports multiple input formats (.dat, .dal, .ou).
- Interactive command-line interface for manual resonance refinement.
- Fast processing: Handles large datasets in seconds.
- Outputs
  - Stabilization diagram plots,
  - Global DOS vs. energy panorama plots,
  - Summary table of detected resonances and their fitted parameters,
  - DOS peak & Lorentzian fit plots.


## How to cite
If you use **DOSmax** in your research, please cite the following article (currently in preparation):

> Johanna Langner, Anjan Sadhukhan, Henryk A. Witek, and Jayanta K. Saha.  
> *An efficient algorithm to determine the resonance parameters of three-body systems using the stabilization method: A case study for natural parity 1,3L^π (L = 0, 1) doubly-excited resonance states of helium*.  
> Manuscript in preparation, 2025.  
> [https://github.com/giogina/DOSmax](https://github.com/giogina/DOSmax)
 
Download citation files here:
- [BibTeX](./cite/DOSmax.bib)
- [RIS](./cite/DOSmax.ris)

## Table of Contents

- [Installation](#installation)
- [Input File Formats](#supported-input-file-formats)
- [Program Workflow](#program-workflow)
  - [Run DOSmax](#run-dosmax)
  - [Stabilization Diagram and Threshold Inputs](#stabilization-diagram-and-threshold-inputs)
  - [DOS Peak Fitting and Resonance Detection](#dos-peak-fitting-and-resonance-detection)
  - [Manual Resonance Refinement](#manual-resonance-refinement)
  - [Output Files](#output-files)
- [License](#license)

# Installation


### System Requirements:
- Python version: 3.10+
- Supported OS: Linux, macOS, and Windows

Follow these steps to install and run **DOSmax** on your system.


## 1. Prerequisites

Before running **DOSmax**, you need:

- **Python 3.10+**
- **`pip`** (Python package installer)
- **`venv`** (Python virtual environment module)
- **`git`** (to clone the repository)

Instructions vary by operating system:

---

### Windows

1. Download and install Python from [https://www.python.org/downloads](https://www.python.org/downloads).
During installation, check:
     - **Add Python to PATH**
     - **Install `pip`**
     - **Install `venv`**

2. Install Git from [https://git-scm.com](https://git-scm.com)

---

### Linux (Debian/Ubuntu)

```bash
    sudo apt update
    sudo apt install python3 python3-pip python3-venv git
```

---
### macOS

1. Install **Homebrew** (if not already available):

 ```bash
     /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
 ```

2. Install Python and git:

 ```bash
    brew install python git
 ```

> Homebrew's Python includes `pip` and `venv` by default.


## 2. Clone the Repository
Clone the **DOSmax** repository from GitHub and navigate into the project directory:

```bash
    git clone https://github.com/giogina/DOSmax.git
```

This downloads the latest version of **DOSmax** locally.

On first startup, if any dependencies are missing, **DOSmax** automatically sets up a virtual environment and downloads the required libraries using `pip`.

Alternatively, if you'd rather install packages system-wide:
```bash
    cd DOSmax/
    pip install -r requirements.txt
```
This will cause **DOSmax** to run using the system-wide Python.



## 3. Run DOSmax

### Windows

In Windows systems, run **DOSmax** using the console command:
```bash
   py path/to/DOSmax.py -f path/to/input_file.dat
```

If that doesn't work, try:
```bash
   python path/to/DOSmax.py -f path/to/input_file.dat
```

### Linux / macOS

You can now execute **DOSmax** with your input file:
```bash
   python3 path/to/DOSmax.py -f path/to/input_file.dat
```
Replace `path/to/input_file.dat` with the actual path to your input file.

On Linux or macOS systems, you can run the script directly (without calling python3) via 
```bash
    path/to/DOSmax.py -f path/to/input_file.dat
```
For this to work, the `./DOSmax.py` script must be executable, which should be the case when cloned from GitHub. If it is not, run:
```bash
   chmod +x DOSmax.py
```

Optionally, to run **DOSmax** from any directory, add it to your system PATH by appending the following line to your `~/.bashrc`, `~/.bash_profile`, or `~/.zshrc`:
```bash
    export PATH="$PATH:/path/to/DOSmax"
```
Then, reload the shell configuration:
```bash
    source ~/.bashrc   # or ~/.bash_profile or ~/.zshrc, depending on your shell
```



# Supported Input File Formats

**DOSmax** accepts input files containing the **diagonalized eigenroot spectrum** for a range of values of the **basis set parameter** $\gamma$. The parser supports three file formats, which are identified automatically by their **file extensions**:


### Tabular Format (`.dat`)
- A **tab-delimited** or **space-delimited** file where:  
  - The **first column** contains **γ values**.  
  - Each **subsequent column** contains **energy values** \( E(\gamma, \text{root}) \) for a specific root.

```text
gamma_1    E_1_root1    E_1_root2    ...    E_1_rootN
gamma_2    E_2_root1    E_2_root2    ...    E_2_rootN
...
gamma_M    E_M_root1    E_M_root2    ...    E_M_rootN
```

- **Notes:**  
  - All rows must have the **same number of columns**.  
  - **Missing values** are **not allowed**—ensure each γ value has corresponding energy entries for all roots.


### Block-Structured Format (`.dal`)

- A file divided into **blocks**, where
  - each block starts with a **single γ value**,
  - followed by the corresponding **energy values** for each root.
  - **Blocks are separated by a blank line** (double newline).  

```text
gamma_1
E_1_root1
E_1_root2
...
E_1_rootN

gamma_2
E_2_root1
E_2_root2
...
E_2_rootN
```


**Notes:**  
- **No extra text** or comments are allowed between blocks.


### Array-Based Format (`.ou`)

- A file containing γ and energy values as arrays, specifically:
  - A **single line** listing all **γ values**.  
  - Subsequent sections, each corresponding to a **root**, listing all associated **energy values** for that root.

```text
gamma_1    gamma_2    ...    gamma_M

# Energies for Root 1
E_1_root1  E_2_root1  ...    E_M_root1

# Energies for Root 2
E_1_root2  E_2_root2  ...    E_M_root2
```

- **Notes:**  
  - **γ values** must be listed **first** and **on one line**, separated by spaces or tabs.  
  - Each **root’s energy array** must be on **a separate line**, following the γ array.
  - The **first energy line** corresponds to **root 1**, the **second** to **root 2**, and so forth.

---

### Input Format Detection
**DOSmax** detects the format **automatically** based on the **file extension**:
- `.dat` → **Tabular format**  
- `.dal` → **Block-structured format**  
- `.ou`  → **Array-based format**  

If an unsupported extension is provided, **DOSmax** will raise an error and list the accepted formats.

### Troubleshooting Input Files:
- **Error:** `ValueError: could not convert string to float`  
  - **Cause:** Non-numeric text or inconsistent formatting.  
  - **Fix:** Check for **hidden characters** or **inconsistent delimiters**.

- **Error:** `IndexError: list index out of range`  
  - **Cause:** Missing energy entries for some roots.  
  - **Fix:** Ensure **complete energy data** for each γ value and root.


# Program workflow

Follow these instructions to run **DOSmax** and interpret its outputs. An example input file, [he_1Po_InfMass.dal](example/he_1Po_InfMass.dal), which was used to create the outputs shown below, is available in the example/ directory of this repository.  

## Run DOSmax

Execute **DOSmax** with your input file:
```bash
python DOSmax.py -f path/to/input_file.dat
```

Replace `path/to/input_file.dat` with the actual path to your input file.

All output files will be written to the directory `path/to/input_file/`, 
therefore it is recommended to first move the input file to a convenient location.


## Stabilization diagram and threshold inputs

After loading and parsing a valid input file, a **stabilization diagram** is displayed:

![Available commands](manual/1_stabilization_diagram.png)

Then, **DOSmax** enters the first interactive stage, which allows users to:
- Visualize the stabilization diagram across various energy ranges.
- Input a list of ionization threshold values, or a nuclear charge implicitly defining these thresholds. 
- Initiate DOS calculations and automatic resonance detection, optionally restricted to a given range.

![Available commands](manual/2_threshold_input.png)

The panorama plot shows log10(DOS) over the entire energy range, with each root colored separately:

![Available commands](manual/stab_he_1Se_singlet_combined_dos_panorama.png)

Selecting the **`r`** command triggers the automatic detection or resonances, which takes about 2 seconds.

## Output Files and Visualization

Upon completion, **DOSmax** generates:

- **results.txt:**  
  - Tab-separated summary of all detected resonances, including:  
    - **Energy** (\(E_r\))  
    - **Width** (\(\Gamma\))  
    - **Amplitude (A)**  
    - **Baseline offset (\(y_0\))**  
    - **Basis set parameter (\(\gamma(E_r)\))**

- **Plot Outputs:**  
  - **Stabilization diagrams**  
  - **Resonance overview plots**  
  - **Lorentzian fit plots** for each detected resonance


## DOS Peak Fitting and Resonance Detection
## Manual Resonance Refinement
## Output Files
# Licence