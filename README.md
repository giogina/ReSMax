# DOSmax
Efficient, automated resonance detection using the Stabilization method.

## Table of Contents

- [Installation](#installation)
- [Usage Guide](#usage-guide)
- [License](#license)

## Features
- Automated resonance detection with minimal user input.
- Supports multiple input formats (.dat, .dal, .ou).
- Interactive command-line interface for manual resonance refinement.
- Fast processing: Handles large datasets in seconds.
- Outputs include stabilization diagrams, DOS plots, and resonance summaries.

# Installation

Follow these steps to install and run **DOSmax** on your system.

---

### 1. Clone the Repository  
Clone the **DOSmax** repository from GitHub and navigate into the project directory:

```bash
git clone https://github.com/giogina/DOSmax.git
cd DOSmax
```

This downloads the latest version of **DOSmax** locally.

---

### 2. (Optional) Install Requirements  
**DOSmax** automatically installs required libraries during the first run. However, you can install them manually:
```bash
pip install -r requirements.txt
```

---

### 3. Run DOSmax  
You can now execute **DOSmax** with your input file:
```bash
python path/to/DOSmax.py -f path/to/input_file.dat
```
Replace `path/to/input_file.dat` with the actual path to your input file.

---

### 4. (Optional) Add DOSmax to Your PATH  
To run **DOSmax** from any directory, add it to your system PATH:

#### On macOS/Linux  
Append the following line to your ````/.bashrc`, ````/.bash_profile`, or ````/.zshrc`:
```bash
export PATH="$PATH:/path/to/DOSmax"
```
Then, reload the shell configuration:
```bash
source ```/.bashrc   # or ```/.zshrc, depending on your shell
```

#### On Windows  
1. Search for **"Environment Variables"** in the Start menu.  
2. Select **"Edit the system environment variables"**.  
3. Click **"Environment Variables"** → Under **System variables**, select **Path** → **Edit**.  
4. Click **"New"** and add the path to the **DOSmax** folder (e.g., `C:\Users\YourName\DOSmax`).  

After setting the PATH, you can simply run:
```bash
DOSmax.py path/to/input_file.dat
```

---

### System Requirements:
- Python version: 3.10  
- Supported OS: Linux, macOS, and Windows  
- No additional installations required: All dependencies are handled automatically by **DOSmax**.


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

### Input Format Selection
**DOSmax** detects the format **automatically** based on the **file extension**:
- `.dat` → **Tabular format**  
- `.dal` → **Block-structured format**  
- `.ou`  → **Array-based format**  

If an unsupported extension is provided, **DOSmax** will raise an error and list the **accepted formats**.


### Important Notes:
- Ensure **no extra blank lines** at the end of files.  
- **File encoding:** UTF-8 is recommended.  
- **Delimiter consistency:**  
  - `.dat` files: Use either **tabs** or **consistent spaces**.  
  - `.dal` and `.ou` files: Use **single spaces** or **tabs** between values.


### Troubleshooting Input Files:
- **Error:** `ValueError: could not convert string to float`  
  - **Cause:** Non-numeric text or inconsistent formatting.  
  - **Fix:** Check for **hidden characters** or **inconsistent delimiters**.

- **Error:** `IndexError: list index out of range`  
  - **Cause:** Missing energy entries for some roots.  
  - **Fix:** Ensure **complete energy data** for each γ value and root.


# Program workflow

Follow these instructions to run **DOSmax** and interpret its outputs.

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
