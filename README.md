# TaxonTableTools2 (TTT)

## Overview

TaxonTableTools2 (TTT) is an easy-to-use graphical software designed for
the analysis and visualization of DNA metabarcoding data. It enables
biologists and researchers without bioinformatics experience to explore
taxonomic datasets quickly, reproducibly, and interactively through a
modern graphical user interface.

TTT focuses on:
* Intuitive biodiversity data exploration
* Reproducible analysis workflows
* Rapid visualization of metabarcoding results
* Accessibility for non-programmers

## Version 2 – Streamlit-based GUI

TaxonTableTools2 represents a complete redesign of the original
TaxonTableTools software.

Because Version 2 was rewritten from the ground up, early releases may
still contain bugs or incomplete features. Stability and functionality
will continue to improve throughout 2026.

Bug reports, feature requests, and suggestions are highly welcome.

## Requirements
* Python3.10 or higher
* Miniconda (recommended)
* Windows (W10 or W11)
* macOS (M1 or M2)
* Linux requires a manual installation via pip

## Conda Installation (recommended)

1.  Download and install [Miniconda](https://www.anaconda.com/download/success)

2.  Open a Miniconda Terminal

Windows: Search for ‘Anaconda Powershell Prompt (Miniconda3)’

macOS: Open a new Terminal window. The prompt should display the (base)
environment.

3.  Download the TTT Environment File

* [Windows](https://github.com/TillMacher/TaxonTableTools2/blob/main/environments/taxontabletools2_env_windows_aarch64.yml)
* [macOS](https://github.com/TillMacher/TaxonTableTools2/blob/main/environments/taxontabletools2_env_macos_aarch64.yml)

4.  Create the TTT Environment (remember to adjust the PATH to match your local file!)
   ```sh
   conda env create -f taxontabletools2_env_windows_aarch64.yml
   ```

5.  Activate the Environment
   ```sh
   conda activate TTT
   ```

6.  Start TaxonTableTools2
   ```sh
   taxontabletools2
   ```

Your browser will automatically open the graphical user interface.

## Manual pip Installation

1. Install python 3.10 or higher (conda environment is recommended)

2. Install taxontabletool2 via pip
   ```sh
   pip install taxontabletools2
   ```

3. Install [scikit-bio](https://scikit.bio/install.html)
   ```sh
   conda install -c conda-forge scikit-bio
   ```
   OR
   ```sh
   pip install scikit-bio
   ```

4. Start TaxonTableTools2
   ```sh
   taxontabletools2
   ```

## Tutorial

Documentation and tutorials are currently under development and will be
released soon.

## Reporting Issues

If you encounter bugs or unexpected behaviour:
* open an issue in the GitHub repository
* contact me directly via email

Community feedback strongly contributes to improving TTT.

## Citation

If you use TaxonTableTools in your research, please cite:

Macher, T. H., Beermann, A. J., & Leese, F. (2021). TaxonTableTools—A
comprehensive, platform-independent graphical user interface software to
explore and visualise DNA metabarcoding data. Molecular Ecology
Resources. https://doi.org/10.1111/1755-0998.13358


## Development Status

TaxonTableTools2 is actively developed and evolving. New analysis
modules, visualization tools, and workflow improvements will be added
continuously.

Contributions and collaborations are welcome.
