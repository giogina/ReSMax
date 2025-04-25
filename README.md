# DOSmax
**Efficient, automated resonance detection using the Stabilization Method.**

**DOSmax** is a Python tool for identifying resonance states in quantum few-body systems, such as helium-like atoms, using the Stabilization Method. This technique analyzes how energy eigenvalues evolve with a systematically varied basis set parameter, revealing metastable bound states as well as resonance states, which appear as characteristic plateaus in the eigenvalue spectrum.

**DOSmax** automates the full workflow: it computes the density of states (DOS) for each root, segments and fits DOS peaks using Lorentzian profiles, groups them into resonances, and extracts resonance parameters (energy and width) from the best-fitting DOS peaks — all with minimal manual input. Results can be interactively inspected and refined via plots and command-line controls, allowing for careful manual validation.

The tool is designed for reproducibility and efficient large-scale analysis of stabilization data, helping researchers extract physical resonance parameters quickly and consistently.

![](example/he_1Po_InfMass/resonance_plots/-0.50000/all_resonances_-0.500.png)

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

<div style="text-align: center;">
  <a href="example/he_1Po_InfMass/resonance_plots/-0.50000/[1]-0.69313909.png" target="_blank">
    <img src="example/he_1Po_InfMass/resonance_plots/-0.50000/[1]-0.69313909.png" alt="Fitted DOS peak" style="max-width: 400px;">
  </a>
</div>



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
- Supported OS: Linux, macOS, and Windows 10+.

Follow these steps to install and run **DOSmax** on your system.


## 1. Prerequisites

Before running **DOSmax**, you need:

- **Python 3.10+**
- **`pip`** (Python package installer)
- **`venv`** (Python virtual environment module)
- **`git`** (to clone the repository)

Instructions vary by operating system:

<details>
<summary><strong>Windows</strong></summary>

1. Download and install Python from [https://www.python.org/downloads](https://www.python.org/downloads).
During installation, check:
     - **Add Python to PATH**
     - **Install `pip`**
     - **Install `venv`**

2. Install Git from [https://git-scm.com](https://git-scm.com)

</details>

<details>
<summary><strong>Linux (Debian/Ubuntu)</strong></summary>

```bash
    sudo apt update
    sudo apt install python3 python3-pip python3-venv git
```

</details>
<details>
<summary><strong>macOS</strong></summary>

1. Install **Homebrew** (if not already available):

 ```bash
     /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
 ```

2. Install Python and git:

 ```bash
    brew install python git
 ```

> Homebrew's Python includes `pip` and `venv` by default.
</details>

## 2. Clone the Repository
Clone the **DOSmax** repository from GitHub and navigate into the project directory:

```bash
    git clone https://github.com/giogina/DOSmax.git
```

This downloads the latest version of **DOSmax** locally.


<details><summary>Optional: Manually install dependencies.</summary>
On first startup, if any dependencies are missing, **DOSmax** automatically sets up a virtual environment and downloads the required libraries using `pip`.

Alternatively, if you would rather install packages system-wide:
```bash
    cd DOSmax/
    pip install -r requirements.txt
```
This will cause **DOSmax** to run using the system-wide Python.
</details>

## 3. Run DOSmax


<details>
<summary><strong>Windows</strong></summary>

In Windows systems, run **DOSmax** using the console command:
```bash
   py path/to/DOSmax.py -f path/to/input_file.dat
```

If that doesn't work, try:
```bash
   python path/to/DOSmax.py -f path/to/input_file.dat
```
</details>

<details>
<summary><strong>Linux / macOS</strong></summary>

You can now execute **DOSmax** with your input file:

```bash
    path/to/DOSmax.py -f path/to/input_file.dat
```
Replace `path/to/input_file.dat` with the actual path to your input file.
For this to work (without calling python3), the `./DOSmax.py` script must be executable, which should be the case when cloned from GitHub. If it is not, first run:
```bash
   chmod +x DOSmax.py
```

<details><summary>Optional: Modify system PATH variable.</summary>

Optionally, to call **DOSmax** from any directory as
```bash
    DOSmax.py -f path/to/input_file.dat
```
add it to your system PATH by appending the following line to your `~/.bashrc`, `~/.bash_profile`, or `~/.zshrc`:
```bash
    export PATH="$PATH:/path/to/DOSmax"
```
Then, reload the shell configuration:
```bash
    source ~/.bashrc   # or ~/.bash_profile or ~/.zshrc, depending on your shell
```

</details>
</details>


# Supported Input File Formats

**DOSmax** accepts input files containing the **diagonalized eigenroot spectrum** for a range of values of the **basis set parameter** $\gamma$. The parser supports three file formats, which are identified automatically by their **file extensions**:

<details>
<summary><strong>Tabular Format (.dat)</strong></summary>

- A **tab-delimited** or **space-delimited** file where:  
  - The **first column** contains **$\gamma$ values**.  
  - Each **subsequent column** contains **energy values** $E(\gamma, \text{root}) $ for a specific root.

```text
gamma_1    E_1_root1    E_1_root2    ...    E_1_rootN
gamma_2    E_2_root1    E_2_root2    ...    E_2_rootN
...
gamma_M    E_M_root1    E_M_root2    ...    E_M_rootN
```
</details>

<details>
<summary><strong>Block-Structured Format (.dal)</strong></summary>

- A file divided into **blocks**, where
  - each block starts with a **single $\gamma$ value**,
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
An example file of this structure, [he_1Po_InfMass.dal](example/he_1Po_InfMass.dal), is available in the example/ directory.
</details>

<details>
<summary><strong>Array-Based Format (.ou)</strong></summary>

- A file containing $\gamma$ and energy values as arrays, specifically:
  - A **single line** listing all **γ values**.  
  - Subsequent sections, each corresponding to a **root**, listing all associated **energy values** for that root.

```text
gamma_1    gamma_2    ...    gamma_M

# Energies for Root 1
E_1_root1  E_2_root1  ...    E_M_root1

# Energies for Root 2
E_1_root2  E_2_root2  ...    E_M_root2
```
</details>

If an unsupported extension is provided, **DOSmax** will raise an error and list the accepted formats.


<details>
<summary><strong>Troubleshooting Input Files</strong></summary>

- **Error:** `ValueError: could not convert string to float`  
  - **Cause:** Non-numeric text or inconsistent formatting.  
  - **Fix:** Check for **hidden characters** or **inconsistent delimiters**.

- **Error:** `IndexError: list index out of range`  
  - **Cause:** Missing energy entries for some roots.  
  - **Fix:** Ensure **complete energy data** for each γ value and root.

</details>

# Program workflow

Follow these instructions to run **DOSmax** and interpret its outputs. In the following, the workflow is demonstrated on the input file [he_1Po_InfMass.dal](example/he_1Po_InfMass.dal), which can be found in the example/ directory of this repository.

All output files are written to a directory of the same name as the input file, in our case `he_1Po_InfMass/`. Therefore, it is recommended to first move the input file to a convenient location before running **DOSmax**.


## Stabilization Diagram and Threshold Inputs

After loading and parsing a valid input file, a **stabilization diagram** is displayed:

![](example/he_1Po_InfMass/stabilization_diagram.png)

Then, **DOSmax** enters the first interactive stage:

![](manual/2_threshold_input.png)

### Thresholds
Here, the ionization threshold of the investigated few-body system can be entered in two different ways: Either by inputting the nuclear charge Z of a helium-like atom, which causes the thresholds to be computed as -Z^2/n^2/2 for n ranging from 1 to 20, or by directly entering a list of thresholds. By default, the ionization thresholds for the helium atom (Z=2) are used.

### DOS Panorama

The panorama plot is a scatter plot of log10(DOS) over the entire energy range for all roots, with each root colored separately:

![](example/he_1Po_InfMass/dos_panorama.png)

This plot gives a clear visual representation of the resonance energies in the form of accumulative peaks, which also makes it clear which peaks are reliable (i.e., clusters of clearly separated higher DOS values), and which are not (in noisy regions below the thresholds).

### Initiate Resonance Detection

Selecting the **`r`** command triggers the automatic detection of resonances, which takes about 2 seconds.

## DOS Peak Fitting and Resonance Detection

### Metastable Bound States
The automatic analysis of the stabilization diagram begins by listing the metastable bound states (MBS), which are detected as local energy minima below the lowest populated ionization threshold. The MBS are listed together with the value of the basis set parameter $\gamma$ at which they appear. If the minimum appears very close to either end of the $\gamma$ range, a warning is issued; these MBS can likely be improved by widening the investigated range of $\gamma$ in the underlying stabilization diagram computation. The printed MBS are also listed in the `results.txt` file.

![](manual/3_MBS.png)

## Manual Resonance Refinement

The DOS curves for each root are now sectioned into individual peaks, which are then fitted to Lorentzian curves, and classified into resonances. The result is first shown graphically in form of a stabilization diagram in which the plateaus corresponding to each DOS peak are highlighted, and colored according to their respective resonance.

![](example/he_1Po_InfMass/resonance_plots/-0.50000/all_resonances_-0.500.png)


## Output Files

Upon completion, **DOSmax** generates and opens a `results.txt` file listing all detected MBS and resonances, grouped by thresholds:

![](manual/5_results_txt.png)

This file lists, for each resonance $r$, the energy $E_r$, Lorentzian fit parameters (Width $\Gamma$, area $A$, offset $y_0$), as well as measures for the fit quality - as sum of squared restudies (SSR), and as relative SSR per point defined by $$\chi_{r} = \frac{\sqrt{SSR}}{N \cdot \rho_{\text{max}}^2}.$$ Additionally, the value of the basis set parameter $\gamma$ at which the chosen DOS peak's maximum occurs is given.

The resonance information is tab-separated, so that it can be copy&pasted into Excel or any other spreadsheet software.

Additionally, plots of the chosen DOS peaks for each resonance, together with their Lorentzian fits, can be generated:
![](manual/6_plot_action.png)

By inputting `p` or `pn`, these plots are generated for all thresholds or the $n^\text{th}$ threshold, respectively. The plot `.png` files are stored in the directory `{file_name}/resonance_plots/{threshold_energy}/`. In our example, the fitted DOS peak for the lowest detected resonance can, after selecting `p1`, be found under `he_1Po_InfMass/resonance_plots/-0.50000/[1]-0.69313909.png`:


<div style="text-align: center;">
  <a href="example/he_1Po_InfMass/resonance_plots/-0.50000/[1]-0.69313909.png" target="_blank">
    <img src="example/he_1Po_InfMass/resonance_plots/-0.50000/[1]-0.69313909.png" alt="Fitted DOS peak" style="max-width: 500px;">
  </a>
</div>

Additionally, for each plotted DOS peak, a `.txt` file of the same name is created, in which the energy and DOS arrays for the DOS data points as well as the fitted Lorentzian curve are given. Using this data, the plot shown above can be re-created in any data visualizer of choice.


# Licence

**DOSmax** is released under the [MIT License](LICENSE).  
You are free to use, modify, and distribute this software with proper attribution.
