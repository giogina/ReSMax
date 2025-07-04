# ReSMax
**Efficient, automated resonance detection using the Stabilization Method.**

**ReSMax** is a Python tool for identifying resonance states in quantum few-body systems, such as helium-like atoms, using the Stabilization Method. This technique analyzes how energy eigenvalues evolve with a systematically varied basis set parameter, revealing metastable bound states as well as resonance states, which appear as characteristic plateaus in the eigenvalue spectrum.

**ReSMax** automates the full workflow: it computes the density of states (DOS) for each root, segments and fits DOS peaks using Lorentzian profiles, groups them into resonances, and extracts resonance parameters (energy and width) from the best-fitting DOS peaks — all with minimal manual input. Results can be interactively inspected and refined via plots and command-line controls, allowing for careful manual validation.

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

<div align="center">
  <a href="https://github.com/giogina/ReSMax/raw/main/example/he_1Po_InfMass/resonance_plots/-0.50000/%5B1%5D-0.69313909.png" target="_blank">
    <img src="https://github.com/giogina/ReSMax/raw/main/example/he_1Po_InfMass/resonance_plots/-0.50000/%5B1%5D-0.69313909.png" alt="Fitted DOS peak" width="600">
  </a>
</div>



## How to cite
If you use **ReSMax** in your research, please cite the following article (currently in preparation):

> Johanna Langner, Anjan Sadhukhan, Henryk A. Witek, and Jayanta K. Saha.  
> *An Algorithm for Automated Extraction of Resonance Parameters from the Stabilization Method*.  
> Manuscript in preparation, 2025.  
> [https://github.com/giogina/ReSMax](https://github.com/giogina/ReSMax)
 
Download citation files here:
- [BibTeX](./cite/ReSMax.bib)
- [RIS](./cite/ReSMax.ris)

## Table of Contents

- [Installation](#installation)
- [Input File Formats](#supported-input-file-formats)
- [Program Workflow](#program-workflow)
  - [Run ReSMax](#run-dosmax)
  - [Stabilization Diagram and Threshold Inputs](#stabilization-diagram-and-threshold-inputs)
  - [DOS Peak Fitting and Resonance Detection](#dos-peak-fitting-and-resonance-detection)
  - [Manual Resonance Refinement](#manual-resonance-refinement)
  - [Output Files](#output-files)
- [License](#license)

# Installation

### System Requirements:
- Python version: 3.10+
- Supported OS: Linux, macOS, and Windows 10+.

Follow these steps to install and run **ReSMax** on your system.


## 1. Prerequisites

Before running **ReSMax**, you need:

- **Python 3.10+**
- **`pip`** (Python package installer)
- **`venv`** (Python virtual environment module)
- **`git`** (to clone the repository)

Instructions vary by operating system:

<details>
<summary><strong>Windows</strong></summary>

When using the [Windows executable](dist/), no preparation is needed.

If you wish to instead run **ReSMax** as a python script:

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

Windows users may download the [Windows executable](dist/).

Otherwise, clone the **ReSMax** repository from GitHub and navigate into the project directory:

```bash
    git clone https://github.com/giogina/ReSMax.git
```

This downloads the latest version of **ReSMax** locally.


<details><summary>Optional: Manually install dependencies.</summary>
On first startup, if any dependencies are missing, **ReSMax** automatically sets up a virtual environment and downloads the required libraries using `pip`.

Alternatively, if you would rather install packages system-wide:
```bash
    cd ReSMax/
    pip install -r requirements.txt
```
This will cause **ReSMax** to run using the system-wide Python.
</details>

## 3. Run ReSMax


<details>
<summary><strong>Windows</strong></summary>

<details>
<summary>Run **ReSMax** from the precompiled Windows executable</summary>

If you downloaded the precompiled `resmax.exe` from GitHub:

1. Open **Command Prompt** (`cmd`) or **PowerShell**
2. Navigate to the folder containing `resmax.exe`
3. Run the program with your input file:

```bash
resmax.exe -f path\to\input_file.dat
```

Replace `path\to\input_file.dat` with the actual path to your input file.  
You can also **drag and drop the file into the terminal** to insert its full path.

> ⚠️ **Note:** Double-clicking `resmax.exe` does nothing unless arguments are provided. Always run it from the command line.

---

<details>
<summary>Optional: Add resmax.exe to your system PATH</summary>

To run `resmax.exe` from any directory without typing its full path:

1. Press `Windows + R`, type `sysdm.cpl`, and press Enter  
2. Go to the **Advanced** tab → Click **Environment Variables**
3. Under **System variables**, select the one called `Path` and click **Edit**
4. Click **New**, then enter the full path to the folder containing `resmax.exe`  
   (e.g., `C:\Users\YourName\Downloads` or wherever you put the file)
5. Click **OK** on all dialogs

Now you can simply run:

```bash
resmax.exe -f input_file.dat
```

from any terminal window.

</details>

</details>



<details>
<summary>Run ReSMax as a Python script</summary>

Run **ReSMax** using the console command:
```bash
   py path/to/resmax.py -f path/to/input_file.dat
```

If that doesn't work, try:
```bash
   python path/to/resmax.py -f path/to/input_file.dat
```

</details>
</details>

<details>
<summary><strong>Linux / macOS</strong></summary>

You can now execute **ReSMax** with your input file:

```bash
    path/to/resmax.py -f path/to/input_file.dat
```
Replace `path/to/input_file.dat` with the actual path to your input file.
For this to work (without calling python3), the `./resmax.py` script must be executable, which should be the case when cloned from GitHub. If it is not, first run:
```bash
   chmod +x resmax.py
```

<details><summary>Optional: Modify system PATH variable.</summary>

Optionally, to call **ReSMax** from any directory as
```bash
    resmax.py -f path/to/input_file.dat
```
add it to your system PATH by appending the following line to your `~/.bashrc`, `~/.bash_profile`, or `~/.zshrc`:
```bash
    export PATH="$PATH:/path/to/resmax"
```
Then, reload the shell configuration:
```bash
    source ~/.bashrc   # or ~/.bash_profile or ~/.zshrc, depending on your shell
```

</details>
</details>


# Supported Input File Formats

**ReSMax** accepts input files containing the **diagonalized eigenroot spectrum** for a range of values of the **basis set parameter** $\gamma$. The parser supports three file formats, which are identified automatically by their **file extensions**:

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

If an unsupported extension is provided, **ReSMax** will raise an error and list the accepted formats.


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

Follow these instructions to run **ReSMax** and interpret its outputs. In the following, the workflow is demonstrated on the input file [he_1Po_InfMass.dal](example/he_1Po_InfMass.dal), which can be found in the example/ directory of this repository.

All output files are written to a directory of the same name as the input file, in our case `he_1Po_InfMass/`. Therefore, it is recommended to first move the input file to a convenient location before running **ReSMax**.


## Stabilization Diagram and Threshold Inputs

After loading and parsing a valid input file, a **stabilization diagram** is displayed:

![](example/he_1Po_InfMass/stabilization_diagram.png)

Then, **ReSMax** enters the first interactive stage:

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

The DOS curves for each root are now automatically sectioned into individual peaks, which are then fitted to Lorentzian curves, and classified into resonances. The result is shown graphically, for each ionization threshold, in form of a stabilization diagram in which the plateaus corresponding to each DOS peak are highlighted, and colored according to their respective resonance.

![](example/he_1Po_InfMass/resonance_plots/-0.50000/all_resonances_-0.500.png)

Each plateau is given a label of the form `iRj`, where `i` specifies the resonance number, and `j` the root number. On the right side of the plot, an accumulative scatter plot diplays $\log_{10}(\text{DOS})$ on the horizontal axis vs. the energy on the vertical axis. This gives additional confirmation as to how well the roots stabilize along each resonance. 

In some cases, the algorithm for automatic selection of best-fitting peaks for each resonance makes a wrong choice. For example, in the plot above, the plateau selected for the resonance number 7 is clearly not a very appropriate one. In such situations, manual refinemens can be made via an interactive loop:

![](manual/4_refinement_menu.png)

To address the issue observed in our example case, we can first plot a zoomed-in version of the overview stabilization diagram by entering:

``plot -0.55 -0.52``

which results in

![](example/he_1Po_InfMass/resonance_plots/-0.50000/resonances_-0.55_-0.52_before.png)

Several observations can be made in this plot: Firstly, some of the plateaus are annotated in red. These are sections in which the energy *descends* with growing $\gamma$, making it impossible to fit the DOS to a Lorentzian peak. This can happen either if the basis set is insufficiently sized to accurately describe this resonance state, or when one encounters fluorescence-active resonance states, as is likely the case e.g. for resonance states 7 and 8 here. **ReSMax** generally chooses fittable DOS peaks to represent resonances, if any are available. Therefore, even though the section annotated by `7R19` can barely be counted as a plateau, it is selected for resonance 7, rather than any of the descending sections observed for roots 25-27. The DOS scatter plot on the right side confirms the presence of very closely spaced resonance states around E=-0.546 and E=-0.527, but the complete absense of any resonance at E=-0.529 (where resonance 7 is currently being placed).

This can be fixed by selecting a different plateau to represent resonance 7, such as the one formed by root 25. To select it, simply enter

``7R25``

This immediately triggers the creation of a new plot with the adjusted peak selection:

![](example/he_1Po_InfMass/resonance_plots/-0.50000/resonances_-0.55_-0.52_after.png)

Note that resonance 7 is now being correctly identified. Its energy, just like the one of resonance 8, is displayed in red on the right side, to signify that a descending section is chosen as its representative. 

We may now verify whether the correct section has been chosen for other resonances, say, resonance 5. This can be achieved by entering

``grid 5``

to obtain a grid plot of DOS vs. energy data points vs. their corresponding Lorentizan fits for all roots forming a plateau at this resonance:


<div align="center">
  <a href="https://github.com/giogina/ReSMax/raw/main/example/he_1Po_InfMass/resonance_plots/-0.50000/%5B5%5D-0.54648347_grid.png" target="_blank">
    <img src="https://github.com/giogina/ReSMax/raw/main/example/he_1Po_InfMass/resonance_plots/-0.50000/%5B5%5D-0.54648347_grid.png" alt="Fitted DOS peak" width="600">
  </a>
</div>

In the top-left panel of this grid, all DOS points are displayed together, to show any drift in the maxima, as well as the relative heights of the DOS values. We can see here that most roots form Lorentzian-like peaks in the same position; jsut the data points from root 25 (cyan) have their maximum values in a different location and fail to form a fittable shape. This is due to the descending nature of the plateau formed by root 25. As all other plateaus correspond to well-behaving DOS peaks, which are fitted well by a Lorentzian curve, we can conclude that the correct one (root 23, green) has been selected.

Closer to each threshold, the resonance states become both more densely packed, and more difficult to describe by any given basis set. This is reflected in the stabilization diagram overview by a more noisy DOS scatter plot on the right sight, where the energy curves on the left side become wavy rather than properly stabilized:


![](example/he_1Po_InfMass/resonance_plots/-0.50000/resonances_-0.53_-0.5-before.png)

Most of the descending sections in the upper part of this plot are spurious, and should not be registered as resonances. We can easily de-select all resonances starting from number 15 by entering

`15+R`

to obtain

![](example/he_1Po_InfMass/resonance_plots/-0.50000/resonances_-0.53_-0.5-after.png)

After being satisfied with all selections for this threshold, we can proceed to the next one by typing

``ok``

Now, just as before, the stabilization diagram overview between the thresholds at E=-0.5 and E=-0.22222 is shown:

![](example/he_1Po_InfMass/resonance_plots/-0.22222/resonances_-0.3_-0.222-before.png)

This area is much noisier, with more wrongly detected resonance states. From the inspection of the right-side scatter plot, it becomes quickly apparent that resonance 48 is entirely spurious; and resonances 51 and up are likely unreliable. We deactivate these resonances by

``48R 51+R``

and obtain

![](example/he_1Po_InfMass/resonance_plots/-0.22222/resonances_-0.3_-0.222-after.png)

This verification process opens a large amount of images; to close all of them, the image viewer process can be terminated by entering `close`. (Note that this only works if the process was started by **ReSMax**.)

The loop over the thresholds can be terminated by typing `end`. This makes sense to do if the presently shown energy range already no longer contains good resonances. No resonances from higher, un-inspected thresholds will be listed in the results.

## Output Files

Upon completion, **ReSMax** generates and opens a `results.txt` file listing all detected MBS and resonances, grouped by thresholds:

![](manual/5_results_txt.png)

This file lists, for each resonance $r$, the energy $E_r$, Lorentzian fit parameters (Width $\Gamma$, area $A$, offset $y_0$), as well as measures for the fit quality - as sum of squared restudies (SSR), and as relative SSR per point defined by $$\chi_{r} = \frac{\sqrt{SSR}}{N \cdot \rho_{\text{max}}^2}.$$ Additionally, the value of the basis set parameter $\gamma$ at which the chosen DOS peak's maximum occurs is given.

The resonance information is tab-separated, so that it can be copy&pasted into Excel or any other spreadsheet software.

Additionally, plots of the chosen DOS peaks for each resonance, together with their Lorentzian fits, can be generated:

![](manual/6_plot_action.png)

By inputting `p` or `pn`, these plots are generated for all thresholds or the $n^\text{th}$ threshold, respectively. The plot `.png` files are stored in the directory `{file_name}/resonance_plots/{threshold_energy}/`. In our example, the fitted DOS peak for the lowest detected resonance can, after selecting `p1`, be found under `he_1Po_InfMass/resonance_plots/-0.50000/[1]-0.69313909.png`:


<div style="text-align: center;">
  <a href="https://github.com/giogina/ReSMax/raw/main/example/he_1Po_InfMass/resonance_plots/-0.50000/%5B1%5D-0.69313909.png" target="_blank">
    <img src="https://github.com/giogina/ReSMax/raw/main/example/he_1Po_InfMass/resonance_plots/-0.50000/%5B1%5D-0.69313909.png" alt="Fitted DOS peak" style="max-width: 400px;">
  </a>
</div>

Additionally, for each plotted DOS peak, a `.txt` file of the same name is created, in which the energy and DOS arrays for the DOS data points as well as the fitted Lorentzian curve are given. Using this data, the plot shown above can be re-created in any data visualizer of choice.


# Licence

**ReSMax** is released under the [MIT License](LICENSE).  
You are free to use, modify, and distribute this software with proper attribution.
