import subprocess
import sys
import os
from pathlib import Path
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
from update_checker import update_check
import importlib.metadata
import threading

## Check for updates
def check_package_update(package_name):
    try:
        # Get the currently installed version
        installed_version = importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        return f"{package_name} is not installed."

    # Check for updates
    res = update_check(package_name, installed_version)
    if res:
        return f"New version of {package_name} is available: {res['latest_version']} (Installed: {installed_version})"
    else:
        return f"{package_name} is up to date (Version: {installed_version})."

# Function to validate the output path
def validate_output_path(path):
    """Check if the given path exists and is a directory."""
    return os.path.isdir(path)

# Function to select a folder
def select_folder():
    global path_to_outdirs
    new_folder = filedialog.askdirectory(title="Select the projects folder")
    if new_folder and validate_output_path(new_folder):
        path_to_outdirs = new_folder
        folder_label.config(text=f"Selected Folder: {path_to_outdirs}")
    else:
        messagebox.showerror("Invalid Folder", "Please select a valid folder.")
        folder_label.config(text=f"Using previous folder: {path_to_outdirs}")

# Function to create a Default_project folder if no projects exist
def ensure_default_project(folder):
    """Check if the folder contains subdirectories (projects) and create Default_project if none exist."""
    sub_folders = [
                    "Alpha_diversity",
                    "Basic_stats",
                    "Beta_diversity",
                    "CCA_plots",
                    "GBIF",
                    "Krona_charts",
                    "Meta_data_table",
                    "NMDS_plots",
                    "Occurrence_analysis",
                    "ParCat_plots",
                    "PCoA_plots",
                    "Per_taxon_statistics",
                    "Perlodes",
                    "Rarefaction_curves",
                    "RAW_tables",
                    "Read_proportions_plots",
                    "Replicate_analysis",
                    "Site_occupancy_plots",
                    "Table_comparison",
                    "Taxon_lists",
                    "TaXon_tables",
                    "TaXon_tables_per_sample",
                    "Taxonomic_resolution_plots",
                    "Taxonomic_richness_plots",
                    "Venn_diagrams"
                ]

    folder_path = Path(folder)
    # Check if the folder is empty or contains only hidden files
    if not any(item.is_dir() or item.is_file() for item in folder_path.iterdir() if not item.name.startswith('.')):
        default_project_path = folder_path.joinpath("Default_project")
        default_project_path.mkdir(parents=True, exist_ok=True)
        # Create subdirectories inside the Default_project folder
        for sub_folder in sub_folders:
            sub_folder_path = default_project_path.joinpath(sub_folder)
            sub_folder_path.mkdir(parents=True, exist_ok=True)
        print(f"Default_project folder created at: {default_project_path}")
    else:
        print(f"Projects found in the selected folder. Default_project folder is not needed.")

# Function to create a new project folder
def create_new_project():
    if not path_to_outdirs or not validate_output_path(path_to_outdirs):
        messagebox.showwarning("Warning", "Please select a valid projects folder first.")
        return

    project_name = simpledialog.askstring("Create New Project", "Enter the name of the new project:")
    if not project_name or not project_name.strip():
        messagebox.showerror("Invalid Name", "Project name cannot be empty.")
        return

    project_path = Path(path_to_outdirs).joinpath(project_name.strip())
    if project_path.exists():
        messagebox.showerror("Error", f"A project named '{project_name}' already exists.")
    else:
        project_path.mkdir(parents=True, exist_ok=True)
        sub_folders = [
            "Alpha_diversity",
            "Basic_stats",
            "Beta_diversity",
            "CCA_plots",
            "GBIF",
            "Krona_charts",
            "Meta_data_table",
            "NMDS_plots",
            "Occurrence_analysis",
            "ParCat_plots",
            "PCoA_plots",
            "Per_taxon_statistics",
            "Perlodes",
            "Rarefaction_curves",
            "RAW_tables",
            "Read_proportions_plots",
            "Replicate_analysis",
            "Site_occupancy_plots",
            "Table_comparison",
            "Taxon_lists",
            "TaXon_tables",
            "TaXon_tables_per_sample",
            "Taxonomic_resolution_plots",
            "Taxonomic_richness_plots",
            "Venn_diagrams"
        ]
        for sub_folder in sub_folders:
            project_path.joinpath(sub_folder).mkdir(parents=True, exist_ok=True)
        messagebox.showinfo("Success", f"New project '{project_name}' created successfully.")

# Function to run the Streamlit app
def run_streamlit_app():
    try:
        # Get the directory of the current script (__main__.py)
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Construct the path to taxontabletools_2.0.py relative to __main__.py
        script_path = os.path.join(script_dir, 'taxontabletools_2.0.py')

        print('Press CTRL + C to close TTT!')
        subprocess.call(['streamlit', 'run', script_path, '--theme.base', 'light', '--server.address', 'localhost'])

    except KeyboardInterrupt:
        sys.exit()

# Function to save the selected folder to the preferences file
def save_folder_to_preferences():
    user_preferences_df.loc[user_preferences_df['Variable'] == 'path_to_outdirs', 'Value'] = path_to_outdirs
    user_preferences_df.to_excel(user_preferences_xlsx, index=False)

# Function to close the tkinter app and continue with the script
def start_app():
    global path_to_outdirs

    # Validate the current output path
    if not path_to_outdirs or not validate_output_path(path_to_outdirs):
        messagebox.showwarning("Warning", "Please select a valid projects folder before continuing.")
        return

    # Save the selected folder to user_preferences.xlsx
    save_folder_to_preferences()

    # Ensure a Default_project exists if no other projects are present
    ensure_default_project(path_to_outdirs)

    # Hide window
    root.withdraw()

    # Run the Streamlit app in a new thread to avoid freezing the GUI
    threading.Thread(target=run_streamlit_app).start()

# Tkinter window setup
root = tk.Tk()
root.title("TaxonTableTools Start Window")
root.geometry("400x300")

# Check for APSCALE package update
update_info = check_package_update('taxontabletools')

# Update label
update_label = tk.Label(root, text=update_info, wraplength=300, justify="left")
update_label.pack(pady=10)

# Paths and user preferences setup
path_to_ttt = Path(__file__).resolve().parent
user_preferences_xlsx = path_to_ttt.joinpath('user_preferences.xlsx')

# Read user preferences from the file
user_preferences_df = pd.read_excel(user_preferences_xlsx).fillna('')
path_to_outdirs = user_preferences_df.loc[user_preferences_df['Variable'] == 'path_to_outdirs', 'Value'].values[0]

# Validate the initial path
if not validate_output_path(path_to_outdirs):
    path_to_outdirs = ''

# Display the current folder
folder_label = tk.Label(root, text=f"Current folder: {path_to_outdirs if path_to_outdirs else 'No folder selected'}")
folder_label.pack(pady=10)

# Button to select a new folder
select_button = tk.Button(root, text="Select Projects Folder", command=select_folder)
select_button.pack(pady=10)

# Button to create a new project
new_project_button = tk.Button(root, text="Create New Project", command=create_new_project)
new_project_button.pack(pady=10)

# Start button to continue with the script
start_button = tk.Button(root, text="Start", command=start_app)
start_button.pack(pady=10)

# Run tkinter main loop
root.mainloop()