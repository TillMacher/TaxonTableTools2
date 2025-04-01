# TaxonTableTools2

[![Downloads](https://pepy.tech/badge/taxontabletools2)](https://pepy.tech/project/taxontabletools2)

## Background

TaxonTableTools2 (TTT) aims to provide easy to use tools for biologists and non-bioinformaticians to analyse and visualize their DNA metabarcoding data quickly and reproducible via a graphical user interface.

TaxonTableTools is an evolving software and there will be bugs and issues at few points. If so, please leave the report in the git repository or drop me an email. Furthermore, new content and functions will be gradually added. Suggestions and recommendations for new features are always welcome!

## Version 2

TaxonTableTools2 has transitioned from PySimpleGUI to Streamlit for its graphical user interface! This shift required a complete rewrite of TTT from the ground up, so the current version may be significantly buggier than previous releases. More stable versions are planned for release throughout 2025.

## Requirements

* Miniconda

## Miniconda Installation

1. Install Miniconda by following the instructions.

2. Open a new Anaconda (Miniconda3) terminal.
   - **Windows**: Type 'Anaconda' in your search bar and select 'Anaconda Powershell Prompt (miniconda3)'.
   - **MacOS**: Open a new terminal. You will see the (base) environment before your user name.

3. Download the respective environment installation file for Windows or [MacOS](https://github.com/TillMacher/TaxonTableTools2/blob/main/environments/taxontabletools_env_macos_aarch64.yml).

4. Install the metabarcoding environment by typing:
   ```sh
   conda env create -f /Users/tillmacher/Downloads/taxontabletools2_env_windows_aarch64.yml
   
5. This should automatically install all dependencies. After the installation, activate the environment:
   ```sh
   conda activate TTT

6. To update TaxonTableTools2 type:
      ```sh
   pip install --upgrade taxontabletools2

## Tutorial

* coming soon

## How to cite

If you use TTT:
* Macher, T.*H., Beermann, A. J., & Leese, F. (2021). TaxonTableTools—A comprehensive, platform*independent graphical user interface software to explore and visualise DNA metabarcoding data. Molecular Ecology Resources. doi: https://doi.org/10.1111/1755*0998.13358

If you create Krona charts, please also cite:
* Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12(1), 385. doi: 10.1186/1471*2105*12*385
