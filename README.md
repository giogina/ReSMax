# DOSmax
Efficient, automated resonance detection using the Stabilization method.

## Features
- Automated resonance detection with minimal user input.
- Supports multiple input formats (.dat, .dal, .ou).
- Interactive command-line interface for manual resonance refinement.
- Fast processing: Handles large datasets in seconds.
- Outputs include stabilization diagrams, DOS plots, and resonance summaries.

## Installation

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
Append the following line to your `~/.bashrc`, `~/.bash_profile`, or `~/.zshrc`:
```bash
export PATH="$PATH:/path/to/DOSmax"
```
Then, reload the shell configuration:
```bash
source ~/.bashrc   # or ~/.zshrc, depending on your shell
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



