import streamlit as st
import pandas as pd
import numpy as np
import subprocess
import glob, sys, time, statistics, os.path, random
from pathlib import Path
import psutil
import plotly.graph_objects as go
import plotly.colors as pc
import streamlit as st
import pandas as pd
import os
import math
import sys
import subprocess
import webbrowser
import glob
import requests, json
from pathlib import Path
import plotly.express as px
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn2_unweighted
from matplotlib_venn import venn3
from matplotlib_venn import venn3_unweighted
from scipy.spatial import distance
from scipy.spatial import ConvexHull
import platform
from skbio.stats.ordination import pcoa
from sklearn.manifold import MDS
from skbio import DistanceMatrix
from skbio.stats.distance import anosim
from pycirclize import Circos
import plotly.io as pio
import plotly.graph_objects as go
from update_checker import UpdateChecker
import importlib.metadata
import statsmodels.api as sm
from scipy import stats
from skbio.diversity.alpha import heip_e
from skbio.diversity.alpha import pielou_e
from skbio.diversity.alpha import shannon
from stqdm import stqdm
from pygbif import occurrences as occ
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from shapely.geometry.polygon import orient
from itertools import product
import shapely.wkt as wkt
from pygbif import maps

########################################################################################################################
## session state cheat sheet

# active_project_path
# table_display_name <- full path to table
# df_taxon_table
# df_taxon_table_uniq <- merge unique taxa to a single entry
# df_metadata_table
# df_samples
# df_traits_table
# df_hashes
# taxa_cols
# taxa_cols2

# TTT_input_path_to_projects
# TTT_height
# TTT_width
# TTT_showlegend
# TTT_template
# TTT_fontsize
# TTT_hash
# TTT_scattersize
# TTT_color1
# TTT_color2
# TTT_colorsequence
# TTT_colorscale

if 'test' == 'pass':
    table_display_name = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test/TaXon_tables/River_mulde_fish_eDNA_2019_merged_NCsub_fish.xlsx')
    df_taxon_table = pd.read_excel(table_display_name).fillna('')
    seq_loc = df_taxon_table.columns.tolist().index('Seq')
    hash_loc = df_taxon_table.columns.tolist().index('Hash')
    species_loc = df_taxon_table.columns.tolist().index('Species')
    df_taxon_cols = df_taxon_table.columns.tolist()[hash_loc+1:species_loc+1]
    trait_cols = ['Hash'] + df_taxon_table.columns.tolist()[species_loc+1:seq_loc]
    df_traits_table = df_taxon_table[trait_cols]
    taxon = df_taxon_table[df_taxon_cols]
    taxon_string = (df_taxon_table[df_taxon_cols].astype(str).agg(";".join, axis=1))
    similarity_loc = df_taxon_table.columns.get_loc("Similarity") + 1
    df_metadata_table = pd.read_excel(table_display_name, sheet_name='Metadata Table').fillna('')
    df_samples = df_metadata_table['Samples'].values.tolist()
    selected_metadata = 'Season'

def default_data_handling(): # just for copying
    # Handle Metadata
    if fig_level != 'Samples':
        df_taxon_table, df_samples = concatenate_by_metadata(fig_level)

    if fig_taxon != 'Hash':
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)
    else:
        df_taxon_table_simple = df_taxon_table.copy()

    if fig_nan == 'Exclude':
        df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']

########################################################################################################################
# Page Setup
st.set_page_config(
    page_title="TaxonTableTools2",
    layout="wide",
    initial_sidebar_state="expanded"
)

########################################################################################################################
## TTT template

TaxonTableTools = go.layout.Template(pio.templates["simple_white"])

# enable grid lines everywhere
TaxonTableTools.layout.xaxis.update(
    showgrid=True,
    gridcolor="lightgrey",
    gridwidth=1
)

TaxonTableTools.layout.yaxis.update(
    showgrid=True,
    gridcolor="lightgrey",
    gridwidth=1
)

# also apply to subplot axes automatically
TaxonTableTools.layout.update(
    xaxis2=dict(showgrid=True, gridcolor="lightgrey"),
    yaxis2=dict(showgrid=True, gridcolor="lightgrey")
)

# register template
pio.templates["TaxonTableTools"] = TaxonTableTools

########################################################################################################################
# help text
help_alphadiv = """
                ### Evenness
                **Evenness** describes how evenly reads (or individuals) are distributed among taxa in a sample.
                
                **High evenness (≈1)**  
                - Taxa have similar abundances  
                - The community is balanced
                
                **Low evenness (≈0)**  
                - A few taxa dominate  
                - Many taxa occur at low abundance
                
                **Interpretation guide**
                - **≈1.0** → very even community  
                - **0.5–0.8** → moderately even  
                - **<0.5** → strong dominance by few taxa  
                
                Evenness measures the **distribution of abundances**, **not the number of taxa**.
                
                ---
                
                ### Shannon Index
                The **Shannon Index** combines **richness (number of taxa)** and **evenness (distribution of abundances)** into a single diversity measure.
                
                **Higher values**
                - Many taxa present
                - Abundances are relatively even
                
                **Lower values**
                - Few taxa present **or**
                - Strong dominance by a few taxa
                
                **Interpretation tips**
                - Shannon values increase with both **more taxa** and **more even communities**
                - Communities with the same richness can have different Shannon values if their **evenness differs**
                """

########################################################################################################################
# general functions

# taxon level plurals
tax_level_plurals = {'Kingdom':'Kingdoms', 'Phylum':'Phyla', 'Class':'Classes', 'Order':'Orders', 'Family':'Families',
                     'Genus':'Genera', 'Species':'Species', 'Hash':'Hashes'}
st.session_state['tax_level_plurals'] = tax_level_plurals
showlabels_options = ['Hide', 'top left', 'top center', 'top right', 'middle left', 'middle center', 'middle right', 'bottom left', 'bottom center', 'bottom right']

def check_package_update(packages):
    for pkg in packages:
        installed_version = importlib.metadata.version(pkg)
        checker = UpdateChecker()
        result = checker.check(pkg, installed_version)
        if result:
            st.sidebar.info(result)
            st.sidebar.info(f'$ pip install --upgrade {pkg}')

def get_package_versions(pkg):
    try:
        version = importlib.metadata.version(pkg)
        return version
    except importlib.metadata.PackageNotFoundError:
        return

def open_folder(folder_path):
    # Get the current operating system
    current_os = platform.system()

    # Open the folder based on the OS
    try:
        if current_os == "Windows":
            subprocess.Popen(f'explorer "{folder_path}"')
        elif current_os == "Darwin":  # macOS
            subprocess.Popen(['open', folder_path])
        else:  # Linux
            subprocess.Popen(['xdg-open', folder_path])
    except Exception as e:
        print(f"Failed to open folder: {e}")

def open_file(path: Path):
    path = str(path)
    if platform.system() == "Darwin":
        subprocess.run(["open", path])
    elif platform.system() == "Windows":
        os.startfile(path)
    else:
        subprocess.run(["xdg-open", path])

def export_taxon_table(table_display_name, df_taxon_table, df_metadata_table, suffix):

    ## collect remaining samples and filter metadata_df
    seq_loc = df_taxon_table.columns.tolist().index('Seq')
    samples = df_taxon_table.columns.tolist()[seq_loc+1::]

    ## Calculate the number of empty OTUs
    empty_otus_count = (~(df_taxon_table[samples] != 0).any(axis=1)).sum()

    ## Remove empty OTUs from dataframe
    mask = (df_taxon_table[samples] != 0).any(axis=1)
    df_taxon_table = df_taxon_table[mask]

    if empty_otus_count != 0:
        st.warning(f'Warning: Removed {empty_otus_count} IDs where all samples were empty (zero reads).')

    ## create a new file name
    if suffix == '':
        output_xlsx = table_display_name
    else:
        output_xlsx = Path(str(Path(table_display_name)).replace('.xlsx', f'_{suffix}.xlsx'))

    ## sort taxon table
    df_taxon_table = df_taxon_table.sort_values(by=['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity'])

    ## export the dataframe
    with pd.ExcelWriter(output_xlsx) as writer:
        df_taxon_table.to_excel(writer, sheet_name='Taxon Table', index=False)
        df_metadata_table.to_excel(writer, sheet_name='Metadata Table', index=False)

def export_plot(folder_name, suffix, fig, lib):

    table_display_name = st.session_state['table_display_name']
    active_project_path = st.session_state['active_project_path']
    os.makedirs(active_project_path / folder_name, exist_ok=True)

    if lib == 'matplot':
        file_pdf = active_project_path / folder_name / f'{suffix}.pdf'
        fig.savefig(file_pdf, dpi=300, bbox_inches="tight")
        st.success(f'Saved plot as pdf!')

    if lib == 'plotly':
        file_pdf = active_project_path / folder_name / f'{suffix}.pdf'
        file_html = active_project_path / folder_name / f'{suffix}.html'
        fig.update_layout(height=st.session_state['TTT_height'], width=st.session_state['TTT_width'])
        fig.write_image(file_pdf)
        fig.write_html(file_html)
        st.success(f'Saved plot as pdf and html!')

def export_table(folder_name, suffix, df, format, index=False):
    table_display_name = st.session_state['table_display_name']
    active_project_path = st.session_state['active_project_path']
    os.makedirs(active_project_path / folder_name, exist_ok=True)

    if format == 'xlsx':
        file_xlsx = active_project_path / folder_name / f'{suffix}.xlsx'
        df.to_excel(file_xlsx, index=index)
        return file_xlsx

    if format == 'csv':
        file_csv = active_project_path / folder_name / f'{suffix}.csv'
        df.to_csv(file_csv, index=index)
        return file_csv

    if format == 'parquet':
        file_parquet = active_project_path / folder_name / f'{suffix}.parquet.snappy'
        df.to_parquet(file_parquet)
        return file_parquet

def rel_taxon_table(df_taxon_table):
    seq_loc = df_taxon_table.columns.tolist().index('Seq')
    df_samples = df_taxon_table.columns.tolist()[seq_loc+1:]
    df = df_taxon_table.copy()
    # Calculate column sums (per sample)
    col_sums = df[df_samples].sum(axis=0)
    # Avoid division by zero
    col_sums = col_sums.replace(0, np.nan)
    # Convert to relative abundance (%)
    df[df_samples] = df[df_samples].div(col_sums, axis=1) * 100
    # Replace NaN (from zero columns) back with 0
    df[df_samples] = df[df_samples].fillna(0)
    return df

def simple_taxon_table(df_taxon_table):

    df = df_taxon_table.copy()

    # Ensure Similarity exists; if not, create placeholder
    if "Similarity" not in df.columns:
        df_taxon_table["Similarity"] = 0

    seq_loc = df.columns.tolist().index('Seq')
    hash_loc = df.columns.tolist().index('Hash')
    species_loc = df.columns.tolist().index('Species')
    df_taxon_cols = df.columns.tolist()[hash_loc+1:species_loc+1]
    df_samples = df.columns.tolist()[seq_loc+1:]

    # Create unique taxon key
    df["taxon_key"] = (df[df_taxon_cols].astype(str).agg("|".join, axis=1))

    def pick_representative(sub):
        # total reads per row
        sub = sub.copy()
        sub["_total_reads"] = sub[df_samples].sum(axis=1)

        # representative row (highest total reads)
        rep_row = sub.loc[sub["_total_reads"].idxmax()]

        return pd.Series(
            {
                "Hash": rep_row["Hash"],
                **{col: rep_row[col] for col in df_taxon_cols},
                "Similarity": sub["Similarity"].max(),
                "N_hashes": len(sub),  # 👈 number of merged hashes
                "Seq": rep_row["Seq"],
                **sub[df_samples].sum().to_dict()
            }
        )

    # Group and collapse
    simple_df = (df.groupby("taxon_key", group_keys=False).apply(pick_representative, include_groups=False).reset_index(drop=True))
    final_cols = (["Hash"] + df_taxon_cols + ["Similarity", "N_hashes", "Seq"] + df_samples)
    simple_df = simple_df[final_cols]

    return simple_df

def get_colors_from_scale(scale_name, n):
    if n == 1:
        return [pc.sample_colorscale(scale_name, [0.5])[0]]

    positions = [i / (n - 1) for i in range(n)]
    return pc.sample_colorscale(scale_name, positions)

def get_colors_from_sequence(sequence_name, n):
    available_colorsequences = st.session_state['available_colorsequences']
    color_sequence = available_colorsequences[sequence_name]

    # repeat the sequence enough times to cover all n
    repeats = (n // len(color_sequence)) + 1
    full_list = color_sequence * repeats
    return full_list[:n]

def colors_to_metadata(df_metadata_table, metadata, colors):
    df_metadata_table = df_metadata_table.copy()
    groups = list(df_metadata_table[metadata].unique())
    colors_dict = {group:color for color, group in zip(colors, groups)}
    res = [colors_dict[sample_color] for sample_color in df_metadata_table[metadata].values.tolist()]
    df_metadata_table['TTT_metadatacolor'] = res
    return df_metadata_table

def concatenate_by_metadata(selected_metadata):

    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()

    # taxonomy columns
    seq_loc = df_taxon_table.columns.get_loc('Seq')
    cols_p1 = df_taxon_table.columns[:seq_loc + 1]

    output_df = df_taxon_table.loc[:, cols_p1].copy()

    # sample → metadata mapping
    sample_to_meta = (
        df_metadata_table
        .loc[df_metadata_table[selected_metadata] != '',
             ['Samples', selected_metadata]]
        .set_index('Samples')[selected_metadata]
    )

    # keep only shared samples
    valid_samples = sample_to_meta.index.intersection(df_taxon_table.columns)

    # column grouping
    grouped_reads = (
        df_taxon_table[valid_samples]
        .T
        .groupby(sample_to_meta.loc[valid_samples])
        .sum()
        .T
    )

    merged_df = pd.concat([output_df, grouped_reads], axis=1)
    updated_samples = merged_df.columns.tolist()[seq_loc+1:]

    return merged_df, updated_samples

def convert_r_p(r, p):
    r = round(r, 3)
    if p >= 0.05:
        p = round(p, 3)
    elif p < 0.05 and p >= 0.001:
        p = round(p, 3)
    elif p < 0.001:
        p = format(p, ".1e")
    return r, p

def check_reads():
    df_taxon_table = st.session_state['df_taxon_table']
    df_samples = st.session_state['df_samples']

    reads = df_taxon_table[df_samples]

    # Convert to numeric, non-numeric values become NaN
    numeric_reads = reads.apply(pd.to_numeric, errors='coerce')

    # Check if any NaN appeared where there was data
    is_numeric = numeric_reads.notna().all().all()

    return is_numeric

########################################################################################################################
# Title
TTT_version = get_package_versions("taxontabletools2")
st.sidebar.markdown(f"# TaxonTableTools v{TTT_version}")
check_package_update(["taxontabletools2"])

########################################################################################################################
# Create Taxon Table
def create_taxon_table():
    # import read table
    selected_read_table = st.session_state['import_folder_files'][st.session_state['selected_read_table']]
    if str(selected_read_table).endswith('.xlsx'):
        read_table_df = pd.read_excel(selected_read_table)
    else:
        read_table_df = pd.read_parquet(selected_read_table)

    # import taxonomy table
    selected_taxonomy_table = st.session_state['import_folder_files'][st.session_state['selected_taxonomy_table']]
    if str(selected_taxonomy_table).endswith('.xlsx'):
        taxonomy_table_df = pd.read_excel(selected_taxonomy_table)
    else:
        taxonomy_table_df = pd.read_parquet(selected_taxonomy_table)

    # verify tables
    # read table
    expected_columns = ['hash', 'sequence']
    for col in expected_columns:
        if col not in read_table_df.columns:
            st.error(f'Could not find the column "{col}" in the Read Table!')
            return

    # taxonomy table
    if st.session_state['taxonomy_table_import_format'] == 'APSCALE':
        expected_columns = ["unique ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Evalue", "Flag", "Ambiguous taxa", "Status"]

        for col in expected_columns:
            if col not in taxonomy_table_df.columns:
                st.error(f'Could not find the column "{col}" in the Taxonomy Table!')
                return

        # check if IDs match
        hashes = sorted(read_table_df['hash'].values.tolist())
        unique_IDs = sorted(taxonomy_table_df['unique ID'].values.tolist())
        if hashes != unique_IDs:
            st.error(f'The hashes between the tables do not match!')
            return

        #  Reorder sequence col
        col = read_table_df.pop("sequence")
        read_table_df.insert(1, "sequence", col)

        # Create "Taxon Table" sheet
        merged_df = (taxonomy_table_df.merge(
            read_table_df,
            left_on="unique ID",
            right_on="hash",
            how="inner").drop(columns=["hash"]))
        merged_df.rename(columns={"unique ID": "Hash"}, inplace=True)
        merged_df.rename(columns={"sequence": "Seq"}, inplace=True)
        merged_df.rename(columns={"evalue": "Evalue"}, inplace=True)

        # Create "Metadata Table" sheet
        seq_loc = merged_df.columns.tolist().index('Seq') + 1
        samples = merged_df.columns.tolist()[seq_loc:]
        metadata_df = pd.DataFrame([[i] + i.split('_') for i in samples])
        metadata_df.rename(columns={0: "Samples"}, inplace=True)

        # Save table to a single Excel file
        to_excel_path = st.session_state['active_project_path'] / 'TaXon_tables' / str(
            st.session_state['taxon_table_creation_name'] + '_taxon_table.xlsx')
        with pd.ExcelWriter(to_excel_path, engine="openpyxl") as writer:
            merged_df.to_excel(writer, sheet_name="Taxon Table", index=False)
            metadata_df.to_excel(writer, sheet_name="Metadata Table", index=False)
        st.success(f'Saved Taxon Table as: {to_excel_path}')

    else:
        expected_columns = ["id", "phylum", "class", "order", "family", "genus", "species", "pct_identity", "status", "records", "records_ratio", "selected_level", "BIN", "flags"]

        for col in expected_columns:
            if col not in taxonomy_table_df.columns:
                st.error(f'Could not find the column "{col}" in the Taxonomy Table!')
                return

        # check if IDs match
        hashes = sorted(read_table_df['hash'].values.tolist())
        unique_IDs = sorted(taxonomy_table_df['id'].values.tolist())
        if hashes != unique_IDs:
            st.error(f'The hashes between the tables do not match!')
            return

        # Create "Taxon Table" sheet
        merged_df = (taxonomy_table_df.merge(
            read_table_df,
            left_on="id",
            right_on="hash",
            how="inner").drop(columns=["hash"]))
        merged_df.rename(columns={"id": "Hash"}, inplace=True)
        merged_df.rename(columns={"sequence": "Seq"}, inplace=True)
        merged_df.rename(columns={"pct_identity": "Similarity"}, inplace=True)
        merged_df.rename(columns={"status": "Status"}, inplace=True)
        merged_df.rename(columns={"records_ratio": "Records ratio"}, inplace=True)
        merged_df.rename(columns={"selected_level": "Selected level"}, inplace=True)
        merged_df.rename(columns={"flags": "Flags"}, inplace=True)
        merged_df.rename(columns={"phylum": "Phylum"}, inplace=True)
        merged_df.rename(columns={"class": "Class"}, inplace=True)
        merged_df.rename(columns={"order": "Order"}, inplace=True)
        merged_df.rename(columns={"family": "Family"}, inplace=True)
        merged_df.rename(columns={"genus": "Genus"}, inplace=True)
        merged_df.rename(columns={"species": "Species"}, inplace=True)
        merged_df.insert(merged_df.columns.get_loc("Hash") + 1, "Kingdom", "")

        # Create "Metadata Table" sheet
        samples = merged_df.columns.tolist()[merged_df.columns.tolist().index('Seq') + 1:]
        metadata_df = pd.DataFrame([[i] + i.split('_') for i in samples])
        metadata_df.rename(columns={0: "Samples"}, inplace=True)

        # Save table to a single Excel file
        to_excel_path = st.session_state['active_project_path'] / 'TaXon_tables' / str(
            st.session_state['taxon_table_creation_name'] + '_taxon_table.xlsx')
        with pd.ExcelWriter(to_excel_path, engine="openpyxl") as writer:
            merged_df.to_excel(writer, sheet_name="Taxon Table", index=False)
            metadata_df.to_excel(writer, sheet_name="Metadata Table", index=False)
        st.success(f'Saved Taxon Table as: {to_excel_path}')

########################################################################################################################
# Tutorial

def TTTutorial():

    st.markdown(
        """
        This short tutorial will guide you through your first steps in **TaxonTableTools2**.
        Follow the steps below to set up and analyze your first project.
        """
    )

    st.markdown("### 1️⃣ Create or Select a Project")

    st.markdown(
        """
        - Enter the path to your **APSCALE projects directory**.
        - TTT2 will automatically create a folder called **`TTT_projects`**.
        - Select your newly created project from the dropdown menu.
        - Click **"Open Active Project"** to access the project folder.

        Each project contains two subfolders:
        - 📂 `Import`
        - 📂 `TaXon_tables`
        """
    )

    st.markdown("### 2️⃣ Import Your Data")

    st.markdown(
        """
        Copy your:
        - **Read Table**
        - **Taxonomy Table**

        into the `Import` folder.

        Supported formats:
        - `.xlsx`
        - `.parquet.snappy`
        """
    )

    path_to_ttt = Path(__file__).resolve().parent
    read_table_df = pd.read_excel(path_to_ttt / 'read_table_example.xlsx').fillna('')
    taxonomy_table_df = pd.read_excel(path_to_ttt / 'taxonomy_table_example.xlsx').fillna('').round(5)

    with st.expander("Example: Read Table"):
        st.dataframe(read_table_df)

    with st.expander("Example: Taxonomy Table"):
        st.dataframe(taxonomy_table_df)

    st.markdown("### 3️⃣ Create a Taxon Table")

    st.markdown(
        """
        Go to **"Create Taxon Table"**:
        - Select your Read Table
        - Select your Taxonomy Table
        - Choose the correct file formats

        TTT2 will automatically merge them into a new **Taxon Table**.

        All Taxon Tables are stored in:
        📂 `TaXon_tables`
        """
    )

    st.markdown("### 4️⃣ Load and Process Tables")

    st.markdown(
        """
        - Load tables using the **"Load Table"** button.
        - If your table does not appear, click **"Refresh"**.

        In the **"Table Processing"** section, you can perform:
        - Replicate merging
        - Negative control subtraction
        - Additional preprocessing steps

        Each processing step creates a new Taxon Table version.
        Remember to refresh and reload before continuing.
        """
    )

    st.markdown("### 5️⃣ Create and Export Plots")

    st.markdown(
        """
        - All visualization modules automatically save plots:
            - as `.pdf`
            - as interactive `.html`

        - The data underlying each plot is also exported as a separate table.
        - Plot appearance (fonts, colors, layout) can be customized via the sidebar.
        """
    )

    st.info(
        "⚠️ TaxonTableTools2 is still under active development. "
        "Please verify your results carefully and report any unexpected behaviour or bugs."
    )

########################################################################################################################
# Basic Stats

def basic_stats_reads():

    """ Calculate the number of reads per sample """
    samples = st.session_state['df_samples']
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    y_values = {i:sum([j for j in df_taxon_table[i].values.tolist() if j != 0]) for i in samples}
    y_values = dict(sorted(y_values.items(), reverse=True, key=lambda item: item[1]))
    x_values = list(y_values.keys())

    max_reads = max(y_values.values())
    min_reads = min(y_values.values())
    avg_reads = round(statistics.mean(y_values.values()), 2)
    try:
        stdev = round(statistics.stdev(y_values.values()), 2)
    except:
        stdev = 0

    fig = go.Figure()
    fig = fig.add_trace(go.Bar(marker_color=st.session_state['TTT_color1'], x=x_values,y=list(y_values.values())))
    fig.update_yaxes(title='Reads', title_font=dict(size=st.session_state['TTT_fontsize']), tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_xaxes(title='Samples', title_font=dict(size=st.session_state['TTT_fontsize']), showticklabels=False, tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_layout(template=st.session_state['TTT_template'], font_size=st.session_state['TTT_fontsize'], title='Distribution of Reads')

    st.plotly_chart(fig, config=st.session_state['TTT_config'])
    st.write(
        f"Max: {max_reads:,.0f} | "
        f"Min: {min_reads:,.0f} | "
        f"Mean: {avg_reads:,.2f} | "
        f"SD: {stdev:,.2f}"
    )

def basic_stats_OTUs():

    """ Calculate the number of OTUs per sample """
    samples = st.session_state['df_samples']
    df_taxon_table = st.session_state['df_taxon_table']
    y_values = {i:len([j for j in df_taxon_table[i].values.tolist() if j != 0]) for i in samples}
    y_values = dict(sorted(y_values.items(), reverse=True, key=lambda item: item[1]))
    x_values = list(y_values.keys())

    max_OTUs, min_OTUs, avg_OTUs = max(y_values.values()), min(y_values.values()), round(statistics.mean(y_values.values()), 2)

    fig = go.Figure()
    fig = fig.add_trace(go.Bar(marker_color=st.session_state['TTT_color1'], x=x_values, y=list(y_values.values())))
    fig.update_layout(template=st.session_state['TTT_template'], font_size=st.session_state['TTT_fontsize'], title=f'Distribution of {st.session_state["TTT_hash"]}')
    fig.update_yaxes(title = st.session_state['TTT_hash'], title_font=dict(size=st.session_state['TTT_fontsize']), tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_xaxes(title='Samples', title_font=dict(size=st.session_state['TTT_fontsize']), showticklabels=False, tickfont=dict(size=st.session_state['TTT_fontsize']))

    st.plotly_chart(fig, config=st.session_state['TTT_config'])
    st.write(
        f"Max: {max_OTUs:,.0f} | "
        f"Min: {min_OTUs:,.0f} | "
        f"Mean: {avg_OTUs:,.2f}"
    )

def top_n_taxa_plot(n=15):
    df_taxon_table = st.session_state['df_taxon_table']
    samples = st.session_state['df_samples']

    # Get top-N taxa
    all_species = [i for i in df_taxon_table['Species'].unique() if i != '']
    species_occurrence_dict = {}
    for species in all_species:
        occurrence = len([i for i in df_taxon_table[df_taxon_table['Species'] == species][samples].sum().values.tolist() if i != 0])
        species_occurrence_dict[species] = occurrence
    species_occurrence_dict = dict(sorted(species_occurrence_dict.items(), key=lambda item: item[1], reverse=True))

    # Plotly bar chart
    x_values = [f'<i>{i}<i>' for i in list(species_occurrence_dict.keys())[:n]]
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=x_values,  # Taxon names on x-axis
        y=list(species_occurrence_dict.values())[:n],
        marker_color=st.session_state['TTT_color1'],
        text=x_values,
        textangle=-90
    ))

    fig.update_layout(template=st.session_state['TTT_template'], font_size=st.session_state['TTT_fontsize'], title=f'Top {n} Species')
    fig.update_yaxes(title = 'Occurrence', title_font=dict(size=st.session_state['TTT_fontsize']), tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_xaxes(title_font=dict(size=st.session_state['TTT_fontsize']), showticklabels=False, title='Species', tickfont=dict(size=st.session_state['TTT_fontsize']*0.6))

    st.plotly_chart(fig, config=st.session_state['TTT_config'])

    st.write(f'{len(all_species)} species across {len(samples)} samples')

def tax_res_plot():
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    taxa_cols = st.session_state['taxa_cols']
    color = st.session_state['TTT_color1']

    y_values = []
    for taxon in taxa_cols:
        n_taxa = len([i for i in df_taxon_table[taxon].unique() if i != ''])
        y_values.append(n_taxa)

    fig = go.Figure()
    fig = fig.add_trace(go.Bar(x=taxa_cols, y=y_values, text=y_values, textposition='outside', cliponaxis=False, marker=dict(color=color)))
    fig.update_layout(template=st.session_state['TTT_template'], font_size=st.session_state['TTT_fontsize'], title=f'Taxonomic resolution')
    fig.update_yaxes(title = "Taxa", title_font=dict(size=st.session_state['TTT_fontsize']), tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_xaxes(title_font=dict(size=st.session_state['TTT_fontsize']), tickfont=dict(size=st.session_state['TTT_fontsize']))

    st.plotly_chart(fig, config=st.session_state['TTT_config'])
    st.write(f'{sum(y_values):,.0f} total taxa')

def collect_sample_stats():
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_samples = st.session_state['df_samples']

    ## information about the table
    n_OTUs = len(df_taxon_table['Hash'])
    n_kingdoms = len(set([i for i in df_taxon_table['Kingdom'].values.tolist() if i != '']))
    n_phyla = len(set([i for i in df_taxon_table['Phylum'].values.tolist() if i != '']))
    n_classes = len(set([i for i in df_taxon_table['Class'].values.tolist() if i != '']))
    n_orders = len(set([i for i in df_taxon_table['Order'].values.tolist() if i != '']))
    n_families = len(set([i for i in df_taxon_table['Family'].values.tolist() if i != '']))
    n_genera = len(set([i for i in df_taxon_table['Genus'].values.tolist() if i != '']))
    n_species = len(set([i for i in df_taxon_table['Species'].values.tolist() if i != '']))
    n_reads = df_taxon_table[df_samples].sum().sum()
    first_row = [n_OTUs, n_kingdoms, n_phyla, n_classes, n_orders, n_families, n_genera, n_species, n_reads]

    samples_info_list = []
    for sample in df_samples:
        sub_df = df_taxon_table[[sample, 'Hash', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
        sub_df = sub_df[sub_df[sample] != 0]

        n_OTUs = len(sub_df['Hash'])
        n_kingdoms = len(set([i for i in sub_df['Kingdom'].values.tolist() if i != '']))
        n_phyla = len(set([i for i in sub_df['Phylum'].values.tolist() if i != '']))
        n_classes = len(set([i for i in sub_df['Class'].values.tolist() if i != '']))
        n_orders = len(set([i for i in sub_df['Order'].values.tolist() if i != '']))
        n_families = len(set([i for i in sub_df['Family'].values.tolist() if i != '']))
        n_genera = len(set([i for i in sub_df['Genus'].values.tolist() if i != '']))
        n_species = len(set([i for i in sub_df['Species'].values.tolist() if i != '']))
        n_reads = sum(sub_df[sample])
        samples_info_list.append([sample, n_OTUs, n_kingdoms, n_phyla, n_classes, n_orders, n_families, n_genera, n_species, n_reads])

    sample_stats_df = pd.DataFrame(samples_info_list,
                                   columns=['Samples', 'Hash', 'Kingdoms', 'Phyla', 'Classes', 'Orders', 'Families', 'Genera', 'Species', 'Reads'])

    # Calculate the mean for each numeric column
    averages = sample_stats_df.select_dtypes(include=[np.number]).mean()
    # Append the averages to the DataFrame using pandas.concat
    sample_stats_df = pd.concat([pd.DataFrame(averages).T, sample_stats_df], ignore_index=True)
    # Set the 'Samples' column of the first row to 'Average'
    sample_stats_df.loc[0, 'Samples'] = 'Average'
    # Replace the index with the 'Samples' column
    sample_stats_df.set_index('Samples', inplace=True)

    # Function to format numbers
    def format_numbers(row):
        if row.name == 'Average':
            return row.round(2).astype(str)
        else:
            return row.astype(int).astype(str)

    # Apply the function to each row
    sample_stats_df = sample_stats_df.apply(format_numbers, axis=1)

    # Insert Dataset information
    sample_stats_df = pd.concat(
        [pd.DataFrame([first_row], columns=sample_stats_df.columns), sample_stats_df],
        ignore_index=False)
    sample_stats_df.index.values[0] = "Dataset"

    st.dataframe(sample_stats_df)

########################################################################################################################
# Table Processing

def sum_if_less_than_X_zeros(row, cutoff):
    if (row != 0).sum() < cutoff:
        return 0
    else:
        return row.sum()

def replicate_merging():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    cutoff_value = st.session_state['cutoff_value']
    missing_mode = st.session_state['missing_replicates_handling']
    suffixes = st.session_state['suffixes']
    df_samples = st.session_state['df_samples']

    # --- derive unique clean sample names ---
    clean_samples = sorted({s.rsplit('_', 1)[0] for s in df_samples})

    # Expected replicate structure
    expected = {
        s: [f"{s}_{suffix}" for suffix in suffixes]
        for s in clean_samples
    }

    # --- detect incomplete replicate groups ---
    incomplete = {
        s: reps for s, reps in expected.items()
        if not all(r in df_taxon_table.columns for r in reps)
    }

    if incomplete:
        st.warning(f"{len(incomplete)} sample(s) have missing replicates.")

        if missing_mode == "Remove":
            cols_to_remove = [
                r for reps in incomplete.values()
                for r in reps if r in df_taxon_table.columns
            ]
            df_taxon_table.drop(columns=cols_to_remove, inplace=True)

    # --- merge complete replicate groups ---
    replicate_cols_to_drop = []

    for s, reps in expected.items():

        # Skip incomplete if we keep them
        if s in incomplete and missing_mode != "Remove":
            continue

        existing_reps = [r for r in reps if r in df_taxon_table.columns]
        if not existing_reps:
            continue

        sub_df = df_taxon_table[existing_reps]

        df_taxon_table[s] = sub_df.apply(
            sum_if_less_than_X_zeros,
            axis=1,
            args=(cutoff_value,)
        )

        replicate_cols_to_drop.extend(existing_reps)

    # Drop replicate columns
    if replicate_cols_to_drop:
        df_taxon_table.drop(columns=replicate_cols_to_drop, inplace=True)

    # --- rebuild metadata based on actual remaining columns ---
    seq_loc = df_taxon_table.columns.get_loc('Seq')
    final_samples = df_taxon_table.columns[seq_loc + 1:]

    df_metadata_table = pd.DataFrame({
        "Samples": final_samples,
        "Placeholder": ""
    })

    # --- export ---
    export_taxon_table(
        table_display_name,
        df_taxon_table,
        df_metadata_table,
        "merged"
    )

    st.success(
        f'Replicates merged → "{table_display_name.name.replace(".xlsx", "_merged.xlsx")}"'
    )

def NC_subtraction():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    df_NC_samples = st.session_state['df_NC_samples']
    NC_handling = st.session_state['NC_handling']
    df_samples = st.session_state['df_samples']
    field_samples = sorted(set(df_samples) - set(df_NC_samples))

    # --- Safety check ---
    if not df_NC_samples:
        st.warning("No negative controls selected.")
        return

    # --- Compute subtraction values per row ---
    if NC_handling == 'Subtract Sum':
        subtract_values = df_taxon_table[df_NC_samples].sum(axis=1)

    elif NC_handling == 'Subtract Max':
        subtract_values = df_taxon_table[df_NC_samples].max(axis=1)

    elif NC_handling == 'Subtract Mean':
        subtract_values = df_taxon_table[df_NC_samples].mean(axis=1)

    elif NC_handling == 'Remove ESVs':
        # Remove rows where any NC has > 0 reads
        mask = df_taxon_table[df_NC_samples].sum(axis=1) == 0
        df_taxon_table = df_taxon_table.loc[mask]
        subtract_values = None

    # Subtract row-wise from field samples
    df_taxon_table[field_samples] = (df_taxon_table[field_samples].sub(subtract_values, axis=0).clip(lower=0))

    # Drop NCs
    df_taxon_table.drop(columns=df_NC_samples, inplace=True, errors="ignore")

    # New metadata df
    df_metadata_table = pd.DataFrame([[i, ''] for i in field_samples], columns=['Samples', 'Placeholder'])

    # Export
    suffix = "NCsub"
    export_taxon_table(
        table_display_name,
        df_taxon_table,
        df_metadata_table,
        suffix
    )

    st.success(f'Negative controls were subtracted and a new table was saved to "{str(table_display_name.name).replace(".xlsx", "_NCsub.xlsx")}"')

def read_based_filter():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_samples = st.session_state['df_samples']

    read_filter = st.session_state['read_filter'] # absolute or relative
    read_filter_mode = st.session_state['read_filter_mode'] # col or row
    read_filter_value = st.session_state['read_filter_value'] # threshold value, can be absolute or relative reads (0-100)

    df = df_taxon_table.copy()
    if read_filter == 'Absolute Reads':
        # Set values below threshold to zero (per sample column)
        df[df_samples] = df[df_samples].where(df[df_samples] >= read_filter_value, 0)
    elif read_filter == 'Relative Reads' and 'column' in read_filter_mode:
        # Convert to relative abundance per column
        col_sums = df[df_samples].sum(axis=0)
        rel_df = df[df_samples].div(col_sums, axis=1) * 100
        df[df_samples] = df[df_samples].where(rel_df >= read_filter_value, 0)
    elif read_filter == 'Relative Reads' and 'row' in read_filter_mode:
        # Relative abundance per row (across all samples)
        total_reads = df[df_samples].sum().sum()
        row_rel = df[df_samples].sum(axis=1) / total_reads * 100
        df = df.loc[row_rel >= read_filter_value]

    # Export
    df_taxon_table = df.copy()
    suffix = "read-filter"
    export_taxon_table(
        table_display_name,
        df_taxon_table,
        st.session_state['df_metadata_table'],
        suffix
    )

    st.success(f'Table was filtered according to read thresholds and a new table was saved to "{str(table_display_name.name).replace(".xlsx", "_readfilter.xlsx")}"')

def read_based_normalisation():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples']
    sub_sample_size = int(st.session_state['sub_sample_size'])

    df = df_taxon_table.copy()
    hash_list = df["Hash"].tolist()

    for sample in df_samples:
        read_df = df[[sample, "Hash"]].copy()
        read_df = read_df[read_df[sample] != 0]
        total_reads = read_df[sample].sum()

        if sub_sample_size <= total_reads:
            # Expand reads into list for subsampling
            read_list = pd.Series(np.repeat(read_df["Hash"].to_list(), read_df[sample].to_list()))
            # Random subsample
            sub_sample = read_list.sample(n=sub_sample_size)
            # Count reads per OTU
            sub_sample_reads = sub_sample.value_counts().to_dict()
            # Reconstruct OTU column in correct order
            df[sample] = [sub_sample_reads.get(otu, 0) for otu in hash_list]
        else:
            # Keep original reads if too few
            df[sample] = df[sample]

    # Export
    suffix = "norm"
    export_taxon_table(
        table_display_name,
        df,
        df_metadata_table,
        suffix
    )

    st.success(f"Table was normalized (subsample size = {sub_sample_size}) and saved to '{str(table_display_name.name).replace('.xlsx', '_norm.xlsx')}'")

def taxonomic_filter():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    taxon_filter_suffix = st.session_state['taxon_filter_suffix']
    taxon_select = st.session_state['taxon_select'] # taxonomic level
    selected_taxa = st.session_state['selected_taxa'] # list of taxa
    taxon_filter_action = st.session_state['taxon_filter_action']

    if len(selected_taxa) == 0:
        st.warning('Please select at least one taxon!')
        return

    # filter table
    if taxon_filter_action == 'Keep': # and remove the rest
        df = df_taxon_table.copy()
        df_filtered = df[df[taxon_select].isin(selected_taxa)]
    elif taxon_filter_action == 'Remove':
        df = df_taxon_table.copy()
        df_filtered = df[~df[taxon_select].isin(selected_taxa)]

    # Export
    export_taxon_table(
        table_display_name,
        df_filtered,
        st.session_state['df_metadata_table'],
        taxon_filter_suffix
    )

    st.success(f'Table was filtered according to the selected taxa and a new table was saved to "{str(table_display_name.name).replace(".xlsx", f"_{taxon_filter_suffix}.xlsx")}"')

def sample_filter():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    sample_filter_suffix = st.session_state['sample_filter_suffix']
    sample_filter_action = st.session_state['sample_filter_action']
    samples_to_filter = st.session_state['samples_to_filter']

    if len(samples_to_filter) == 0:
        st.warning('Please select at least one sample!')
        return

    # filter table
    df = df_taxon_table.copy()
    seq_loc = df.columns.tolist().index('Seq')
    keep_cols = df.columns.tolist()[:seq_loc+1]

    if sample_filter_action == 'Keep': # and remove the rest
        keepers = keep_cols + samples_to_filter
        df_filtered = df[keepers]
    elif sample_filter_action == 'Remove': # and keep the rest
        all_samples = set(st.session_state['df_samples'])
        keepers = all_samples - set(samples_to_filter)
        mask = list(keep_cols) + list(keepers)
        df_filtered = df[mask]

    # Update Metadata Table
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_metadata_table = df_metadata_table[df_metadata_table['Samples'].isin(keepers)]

    # Export
    export_taxon_table(
        table_display_name,
        df_filtered,
        df_metadata_table,
        sample_filter_suffix
    )

    st.success(f'Table was filtered according to the selected samples and a new table was saved to "{str(table_display_name.name).replace(".xlsx", f"_{sample_filter_suffix}.xlsx")}"')

def trait_filter():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    trait_filter_suffix = st.session_state['trait_filter_suffix']
    selected_trait = st.session_state['selected_trait'] # taxonomic level
    selected_trait_values = st.session_state['selected_trait_values'] # list of taxa
    trait_filter_action = st.session_state['trait_filter_action']

    if len(selected_trait) == 0:
        st.warning('Please select at least one taxon!')
        return

    if st.session_state['trait_type'] == 'string':
        # filter table
        if trait_filter_action == 'Keep': # and remove the rest
            df = df_taxon_table.copy()
            df_filtered = df[df[selected_trait].isin(selected_trait_values)]
        elif trait_filter_action == 'Remove':
            df = df_taxon_table.copy()
            df_filtered = df[~df[selected_trait].isin(selected_trait_values)]
    else:
        # filter table
        if trait_filter_action == 'Keep': # and remove the rest
            df = df_taxon_table.copy()
            df_filtered = df[df[selected_trait] >= (selected_trait_values)]
        elif trait_filter_action == 'Remove':
            df = df_taxon_table.copy()
            df_filtered = df[~df[selected_trait] >= (selected_trait_values)]

    # Export
    export_taxon_table(
        table_display_name,
        df_filtered,
        st.session_state['df_metadata_table'],
        trait_filter_suffix
    )

    st.success(f'Table was filtered according to the selected traits and a new table was saved to "{str(table_display_name.name).replace(".xlsx", f"_{trait_filter_suffix}.xlsx")}"')

########################################################################################################################
# Table conversion

def merge_taxon_tables():

    # General
    TTT_hash = st.session_state['TTT_hash']

    # Load second table
    df2_taxon_table_file = st.session_state['df2_taxon_table_file']

    if str(df2_taxon_table_file) == 'None':
        return
    elif str(df2_taxon_table_file.name).endswith('.xlsx'):
        df2_taxon_table = pd.read_excel(df2_taxon_table_file).fillna('')
        df2_metadata_table = pd.read_excel(df2_taxon_table_file, sheet_name='Metadata Table').fillna('')
    else:
        st.warning('Please provide a table in .xlsx format!')
        return

    # Load first table
    df1_taxon_table = st.session_state['df_taxon_table']
    df1_metadata_table = st.session_state['df_metadata_table']
    df1_samples = st.session_state['df_samples']
    table_display_name = st.session_state['table_display_name']
    table_merge_suffix = st.session_state['table_merge_suffix']

    # Identify column structure
    species_loc1 = df1_taxon_table.columns.get_loc('Species')
    seq_loc1 = df1_taxon_table.columns.get_loc('Seq')
    trait_cols1 = list(df1_taxon_table.columns[species_loc1+2:seq_loc1])
    df1_traits = df1_taxon_table[['Hash'] + trait_cols1]
    df1_taxon_table = df1_taxon_table.drop(columns=trait_cols1)

    species_loc2 = df2_taxon_table.columns.get_loc('Species')
    seq_loc2 = df2_taxon_table.columns.get_loc('Seq')
    trait_cols2 = list(df2_taxon_table.columns[species_loc2+2:seq_loc2])
    df2_samples = df2_taxon_table.columns[seq_loc2+1:]
    df2_traits = df2_taxon_table[['Hash'] + trait_cols2]
    df2_taxon_table = df2_taxon_table.drop(columns=trait_cols2)

    ## calculate shared traits
    shared_traits = list(set(trait_cols1) & set(trait_cols2))

    if len(set(df1_samples) & set(df2_samples)) != 0:
        st.error('Error: Cannot merge tables with overlapping sample names!')
    else:
        ## Collect some stats
        ESVs_1 = set(df1_taxon_table['Hash'].values.tolist())
        ESVs_2 = set(df2_taxon_table['Hash'].values.tolist())
        a_only = len(ESVs_1 - ESVs_2)
        shared = len(ESVs_1 & ESVs_2)
        b_only = len(ESVs_2 - ESVs_1)
        total = a_only + shared + b_only
        all_hashes = sorted(set(df1_taxon_table['Hash'].values.tolist() + df2_taxon_table['Hash'].values.tolist()))

        st.markdown(f"""
        **Shared ESVs:** {shared}  
        **Original exclusive {TTT_hash}:** {a_only}  
        **Added exclusive {TTT_hash}:** {b_only}  
        **Total ESVs:** {total}
        """)

        merged_taxon_table_df = pd.DataFrame()
        for hash in stqdm(all_hashes, desc=f'Looping through {TTT_hash}'):
            new_row = pd.DataFrame()

            # Case 1: Present in both datasets: Merge them and select the taxonomy with the higher similarity
            if hash in ESVs_1 and hash in ESVs_2:
                # data from table 1
                row1 = df1_taxon_table.loc[df1_taxon_table['Hash'] == hash]
                data1 = row1.iloc[:, :10]
                similarity1 = row1['Similarity'].values.tolist()[0]
                reads1 = row1[df1_samples]
                traits1 = df1_traits.loc[df1_traits['Hash'] == hash][shared_traits]
                traits1['ESV merge'] = ['Shared']

                # data from table 2
                row2 = df2_taxon_table.loc[df2_taxon_table['Hash'] == hash]
                data2 = row2.iloc[:, :10]
                similarity2 = row2['Similarity'].values.tolist()[0]
                reads2 = row2[df2_samples]
                traits2 = df2_traits.loc[df2_traits['Hash'] == hash][shared_traits]

                if data1.values.tolist() == data2.values.tolist():
                    # Find the index of the column "Seq"
                    col_index = data1.columns.get_loc("Seq")
                    # Insert the new DataFrame columns before "Seq"
                    data1_traits = pd.concat([data1.iloc[:, :col_index], traits1, data1.iloc[:, col_index:]], axis=1)
                    new_row = pd.concat([data1_traits.reset_index(drop=True),
                                         reads1.reset_index(drop=True),
                                         reads2.reset_index(drop=True)], axis=1)

                elif similarity1 > similarity2:
                    # Find the index of the column "Seq"
                    col_index = data1.columns.get_loc("Seq")
                    # Insert the new DataFrame columns before "Seq"
                    data1_traits = pd.concat([data1.iloc[:, :col_index], traits1, data1.iloc[:, col_index:]], axis=1)
                    new_row = pd.concat([data1_traits.reset_index(drop=True),
                                         reads1.reset_index(drop=True),
                                         reads2.reset_index(drop=True)], axis=1)
                else:
                    # Find the index of the column "Seq"
                    col_index = data2.columns.get_loc("Seq")
                    # Insert the new DataFrame columns before "Seq"
                    data2_traits = pd.concat([data2.iloc[:, :col_index], traits2, data2.iloc[:, col_index:]], axis=1)
                    new_row = pd.concat([data2_traits.reset_index(drop=True),
                                         reads1.reset_index(drop=True),
                                         reads2.reset_index(drop=True)], axis=1)

            # Case 2: Present in dataset 1: Fill up the other
            elif hash in ESVs_1:
                # data from table 1
                row1 = df1_taxon_table.loc[df1_taxon_table['Hash'] == hash]
                data1 = row1.iloc[:, :10]
                reads1 = row1[df1_samples]
                traits1 = df1_traits.loc[df1_traits['Hash'] == hash][shared_traits]
                traits1['ESV merge'] = ['Original']

                # data from table 2
                reads2 = pd.DataFrame([[0] * len(df2_samples)], columns=df2_samples)

                # Find the index of the column "Seq"
                col_index = data1.columns.get_loc("Seq")
                # Insert the new DataFrame columns before "Seq"
                data1_traits = pd.concat([data1.iloc[:, :col_index], traits1, data1.iloc[:, col_index:]], axis=1)
                new_row = pd.concat([data1_traits.reset_index(drop=True),
                                     reads1.reset_index(drop=True),
                                     reads2.reset_index(drop=True)], axis=1)

            # Case 2: Present in dataset 1: Fill up the other
            elif hash in ESVs_2:
                # data from table 1
                row2 = df2_taxon_table.loc[df2_taxon_table['Hash'] == hash]
                data2 = row2.iloc[:, :10]
                reads2 = row2[df2_samples]
                traits2 = df2_traits.loc[df2_traits['Hash'] == hash][shared_traits]
                traits2['ESV merge'] = ['Added']

                # data from table 2
                reads1 = pd.DataFrame([[0] * len(df1_samples)], columns=df1_samples)

                # Find the index of the column "Seq"
                col_index = data2.columns.get_loc("Seq")
                # Insert the new DataFrame columns before "Seq"
                data2_traits = pd.concat([data2.iloc[:, :col_index], traits2, data2.iloc[:, col_index:]], axis=1)
                new_row = pd.concat([data2_traits.reset_index(drop=True),
                                     reads1.reset_index(drop=True),
                                     reads2.reset_index(drop=True)], axis=1)

            else:
                print('Whoops that should not happen!')

            ## add to newly merged df
            merged_taxon_table_df = pd.concat([merged_taxon_table_df, new_row], ignore_index=True)

    ## the merged_df is the new Taxon Table
    ## now create the Metadata Table
    ## calculate shared metadata
    shared_metadata = list(set(df1_metadata_table.columns.tolist()) & set(df2_metadata_table.columns.tolist()))
    export_metadata_table_df = pd.concat([df1_metadata_table[shared_metadata], df2_metadata_table[shared_metadata]], ignore_index=True)
    export_metadata_table_df['ESV merge'] = ['Original'] * len(df1_metadata_table) + ['Added'] * len(df2_metadata_table)

    ## export df
    export_taxon_table(table_display_name, merged_taxon_table_df, export_metadata_table_df, table_merge_suffix)

    st.success(f'Tables were successfully merged!')

def simplify_table(save=True):

    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_taxon_cols = st.session_state['df_taxon_cols']
    df_samples = st.session_state['df_samples']  # sample columns
    df_metadata_table = st.session_state['df_metadata_table']

    # Ensure Similarity exists; if not, create placeholder
    if "Similarity" not in df_taxon_table.columns:
        df_taxon_table["Similarity"] = 0

    # Create unique taxon key
    df_taxon_table["taxon_key"] = (df_taxon_table[df_taxon_cols].astype(str).agg("|".join, axis=1))

    def pick_representative(sub):
        # total reads per row
        sub = sub.copy()
        sub["_total_reads"] = sub[df_samples].sum(axis=1)

        # representative row (highest total reads)
        rep_row = sub.loc[sub["_total_reads"].idxmax()]

        return pd.Series(
            {
                "Hash": rep_row["Hash"],
                **{col: rep_row[col] for col in df_taxon_cols},
                "Similarity": sub["Similarity"].max(),
                "N_hashes": len(sub),  # 👈 number of merged hashes
                "Seq": rep_row["Seq"],
                **sub[df_samples].sum().to_dict()
            }
        )

    # Group and collapse
    collapsed_df = (df_taxon_table
                    .groupby("taxon_key", group_keys=False)
                    .apply(pick_representative)
                    .reset_index(drop=True)
                    )

    # Final column order
    final_cols = (
            ["Hash"]
            + df_taxon_cols
            + ["Similarity", "N_hashes", "Seq"]
            + df_samples
    )

    collapsed_df = collapsed_df[final_cols]

    if save == True:
        ## export the dataframe
        simplified_table_xlsx = str(Path(table_display_name)).replace('.xlsx', '_simplified.xlsx')
        with pd.ExcelWriter(simplified_table_xlsx) as writer:
            collapsed_df.to_excel(writer, sheet_name='Taxon Table', index=False)
            df_metadata_table.to_excel(writer, sheet_name='Metadata Table', index=False)
        st.success(f'Table was simplified and saved to "{Path(simplified_table_xlsx).name}"')
    else:
        return collapsed_df

def add_traits_from_file():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    trait_import_df = st.session_state['trait_import_df'].copy()
    trait_import_taxon = st.session_state['trait_import_taxon']
    trait_cols = trait_import_df.columns.tolist()[1:]

    if trait_import_taxon not in st.session_state['taxa_cols']:
        st.error(f'Please choose a vaild taxon as first column: {", ".join(st.session_state["taxa_cols"])}')
        return

    # Merge first
    df_merged = df_taxon_table.merge(
        trait_import_df,
        how='left',
        left_on=trait_import_taxon,
        right_on=trait_import_taxon
    )

    # Reorder columns: all columns up to 'Seq', then traits, then 'Seq', then rest
    cols = df_merged.columns.tolist()
    seq_idx = cols.index('Seq')
    # Exclude traits already in the list
    new_cols = cols[:seq_idx] + trait_cols + ['Seq'] + [c for c in cols[seq_idx + 1:] if c not in trait_cols]
    df_merged = df_merged[new_cols]

    # Export
    export_taxon_table(
        table_display_name,
        df_merged,
        st.session_state['df_metadata_table'],
        ''
    )
    st.session_state['df_taxon_table'] = df_merged.copy()
    st.success(f'Traits were added and the original table was updated!')

def sort_samples():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table']
    sorted_samples = [item[0] for item in df_metadata_table[['Samples']].values.tolist()]

    seq_loc = df_taxon_table.columns.tolist().index('Seq')
    new_col_sort = df_taxon_table.columns.tolist()[:seq_loc+1] + sorted_samples
    df_sorted = df_taxon_table[new_col_sort]

    # Export
    export_taxon_table(
        table_display_name,
        df_sorted,
        st.session_state['df_metadata_table'],
        ''
    )
    st.session_state['df_taxon_table'] = df_sorted.copy()
    st.success(f'Samples were sorted and the original table was updated!')

def sample_rename():
    table_display_name = st.session_state['table_display_name']
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    sample_rename_metadata = st.session_state['sample_rename_metadata']
    rename_dict = {i: j for i, j in df_metadata_table[['Sample', sample_rename_metadata]].values.tolist()}

    df_taxon_table = df_taxon_table.rename(columns=rename_dict)
    df_metadata_table['Old Sample Names'] = list(rename_dict.keys())
    df_metadata_table['Sample'] = list(rename_dict.values())

    # Export
    export_taxon_table(
        table_display_name,
        df_taxon_table,
        df_metadata_table,
        ''
    )
    st.session_state['df_taxon_table'] = df_taxon_table.copy()
    st.session_state['df_metadata_table'] = df_metadata_table.copy()
    st.success(f'Samples were sorted and the original table was updated!')

def export_fasta():
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    hashes = df_taxon_table['Hash'].values.tolist()
    sequences = df_taxon_table['Seq'].values.tolist()

    output_folder = st.session_state['active_project_path'] / 'FASTA'
    os.makedirs(output_folder, exist_ok=True)
    table_display_name = st.session_state['table_display_name']
    output_file = output_folder / f'{table_display_name.stem}_hashes.fasta'

    with open(output_file, 'w') as f:
        for header, seq in zip(hashes, sequences):
            f.write(f'>{header}\n')
            f.write(f'{seq}\n')

    st.success(f'Hashes were saved as .FASTA file!')

########################################################################################################################
# Sample Comparison

def venn_diagram():

    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    venn_taxon_level = st.session_state['venn_taxon_level']
    venn_metadata = st.session_state['venn_metadata']
    venn_display = st.session_state['venn_display']
    venn_weighted = st.session_state['venn_weighted']

    venn_categories = [i for i in sorted(df_metadata_table[venn_metadata].unique()) if i != '']

    if len(venn_categories) == 2:
        name_a = venn_categories[0]
        name_b = venn_categories[1]

        samples_a = df_metadata_table[df_metadata_table[venn_metadata] == name_a]['Samples'].values.tolist()
        if samples_a == str:
            samples_a = [samples_a]
        df_out = df_taxon_table[[venn_taxon_level]].copy()
        df_out["Sum_reads"] = df_taxon_table[samples_a].sum(axis=1)
        species_a = set([i for i in df_out[df_out['Sum_reads'] != 0][venn_taxon_level].unique() if i != ''])

        samples_b = df_metadata_table[df_metadata_table[venn_metadata] == name_b]['Samples'].values.tolist()
        if samples_b == str:
            samples_b = [samples_b]
        df_out = df_taxon_table[[venn_taxon_level]].copy()
        df_out["Sum_reads"] = df_taxon_table[samples_b].sum(axis=1)
        species_b = set([i for i in df_out[df_out['Sum_reads'] != 0][venn_taxon_level].unique() if i != ''])

        a_only = species_a - species_b
        n_a_only = len(a_only)
        shared = species_a & species_b
        n_shared = len(shared)
        b_only = species_b - species_a
        n_b_only = len(b_only)

        fig = plt.figure()

        if venn_weighted == 'Weighted':
            v = venn2(
                subsets=(n_a_only, n_b_only, n_shared),
                set_labels=(name_a, name_b),
                set_colors=(st.session_state['TTT_color1'], st.session_state['TTT_color2']),
                alpha=0.6
            )
        else:
            v = venn2_unweighted(
                subsets=(n_a_only, n_b_only, n_shared),
                set_labels=(name_a, name_b),
                set_colors=(st.session_state['TTT_color1'], st.session_state['TTT_color2']),
                alpha=0.6
            )

        total = n_a_only + n_b_only + n_shared

        # Access labels safely
        label_10 = v.get_label_by_id('10')  # A only
        label_01 = v.get_label_by_id('01')  # B only
        label_11 = v.get_label_by_id('11')  # Shared

        if venn_display == 'Values':
            # Keep absolute values (default)
            pass

        elif venn_display == 'Rel. Values':
            if total > 0:
                if label_10:
                    label_10.set_text(f"{n_a_only / total:.1%}")
                if label_01:
                    label_01.set_text(f"{n_b_only / total:.1%}")
                if label_11:
                    label_11.set_text(f"{n_shared / total:.1%}")

        elif venn_display == 'Both':
            if total > 0:
                if label_10:
                    label_10.set_text(
                        f"{n_a_only}\n({n_a_only / total:.1%})"
                    )
                if label_01:
                    label_01.set_text(
                        f"{n_b_only}\n({n_b_only / total:.1%})"
                    )
                if label_11:
                    label_11.set_text(
                        f"{n_shared}\n({n_shared / total:.1%})"
                    )


        # Font sizes
        for text in v.subset_labels:
            if text:
                text.set_fontsize(st.session_state['TTT_fontsize'])
        for text in v.set_labels:
            if text:
                text.set_fontsize(st.session_state['TTT_fontsize'])
        if venn_taxon_level == 'Hash':
            venn_taxon_level = st.session_state['TTT_hash'][:-1]
        plt.title(f"{venn_taxon_level} Overlap", fontsize=st.session_state['TTT_fontsize'])
        st.pyplot(fig)

        export_plot('Venn_diagrams', f'{name}_{venn_metadata}_{venn_taxon_level}', plt, 'matplot')

    elif len(venn_categories) == 3:
        name_a, name_b, name_c = venn_categories

        def get_species(name):
            samples = df_metadata_table[
                df_metadata_table[venn_metadata] == name
                ]['Samples'].values.tolist()

            df_tmp = df_taxon_table[[venn_taxon_level]].copy()
            df_tmp["Sum_reads"] = df_taxon_table[samples].sum(axis=1)

            return set(
                i for i in
                df_tmp[df_tmp["Sum_reads"] != 0][venn_taxon_level].unique()
                if i != ''
            )

        species_a = get_species(name_a)
        species_b = get_species(name_b)
        species_c = get_species(name_c)

        # Region calculations
        a_only = species_a - species_b - species_c
        b_only = species_b - species_a - species_c
        c_only = species_c - species_a - species_b

        ab = (species_a & species_b) - species_c
        ac = (species_a & species_c) - species_b
        bc = (species_b & species_c) - species_a

        abc = species_a & species_b & species_c

        fig = plt.figure()

        if venn_weighted == 'Weighted':
            v = venn3(
                subsets=(
                    len(a_only),
                    len(b_only),
                    len(ab),
                    len(c_only),
                    len(ac),
                    len(bc),
                    len(abc)
                ),
                set_labels=(name_a, name_b, name_c),
                set_colors=(
                    st.session_state['TTT_color1'],
                    st.session_state['TTT_color2'],
                    st.session_state['TTT_color3']
                ),
                alpha=0.6
            )
        else:
            v = venn3_unweighted(
                subsets=(
                    len(a_only),
                    len(b_only),
                    len(ab),
                    len(c_only),
                    len(ac),
                    len(bc),
                    len(abc)
                ),
                set_labels=(name_a, name_b, name_c),
                set_colors=(
                    st.session_state['TTT_color1'],
                    st.session_state['TTT_color2'],
                    st.session_state['TTT_color3']
                ),
                alpha=0.6
            )

        total = (
                len(a_only) + len(b_only) + len(c_only)
                + len(ab) + len(ac) + len(bc)
                + len(abc)
        )

        # Region IDs for venn3:
        region_map = {
            '100': len(a_only),
            '010': len(b_only),
            '110': len(ab),
            '001': len(c_only),
            '101': len(ac),
            '011': len(bc),
            '111': len(abc)
        }

        if venn_display != 'Values' and total > 0:
            for region_id, value in region_map.items():
                label = v.get_label_by_id(region_id)
                if label:
                    if venn_display == 'Rel. Values':
                        label.set_text(f"{value / total:.1%}")
                    elif venn_display == 'Both':
                        label.set_text(
                            f"{value}\n({value / total:.1%})"
                        )

        # Font size
        for text in v.subset_labels:
            if text:
                text.set_fontsize(st.session_state['TTT_fontsize'])

        for text in v.set_labels:
            if text:
                text.set_fontsize(st.session_state['TTT_fontsize'])

        if venn_taxon_level == 'Hash':
            venn_taxon_level = st.session_state['TTT_hash'][:-1]

        plt.title(
            f"{venn_taxon_level} Overlap",
            fontsize=st.session_state['TTT_fontsize']
        )

        st.pyplot(fig)

        export_plot(
            'Venn_diagrams',
            f'{name}_{venn_metadata}_{venn_taxon_level}',
            plt,
            'matplot'
        )

########################################################################################################################
# Read Proportions

def readdist_diagram():

    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']

    # Tables
    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    # Options
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples'].copy()
    fig_color = st.session_state['readdist_color'] # color sequence or scale
    fig_mode = st.session_state['readdist_mode'] # relative or absolute reads
    fig_level = st.session_state['readdist_level'] # sample or metadata
    fig_taxon = st.session_state['readdist_taxon'] # taxonomic level
    fig_nan = st.session_state['readdist_nan'] # include or exclude

    # Handle Metadata
    if fig_level != 'Samples':
        df_taxon_table, df_samples = concatenate_by_metadata(fig_level)
    df_taxon_table_simple = simple_taxon_table(df_taxon_table)

    # Analysis
    if fig_nan == 'Exclude':
        df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']
    df_taxon_table_simple_s = (df_taxon_table_simple[[fig_taxon] + df_samples].groupby(fig_taxon, as_index=False).sum())
    all_taxa = sorted(df_taxon_table_simple[fig_taxon].unique())[::-1]

    # Output Plot
    # Colors
    if fig_color == 'Color Scale':
        colors = get_colors_from_scale(TTT_colorscale, len(all_taxa))
    if fig_color == 'Color Sequence':
        colors = get_colors_from_sequence(TTT_colorsequence, len(all_taxa))

    fig = go.Figure()
    if fig_mode == 'Relative Reads' or fig_mode == 'Absolute Reads':
        res = []
        for c, taxon in enumerate(all_taxa):
            y_values = df_taxon_table_simple_s[df_taxon_table_simple_s[fig_taxon] == taxon][df_samples].values.tolist()[0]
            x_values = df_samples
            fig.add_trace(go.Bar(x=x_values, y=y_values, name=taxon, marker=dict(color=colors[c])))
            res.append(y_values)
    else:
        res = []
        for c, taxon in enumerate(all_taxa):
            sub_df = df_taxon_table[df_taxon_table[fig_taxon] == taxon][df_samples]
            y_values = [len([i for i in sub_df[sample].values.tolist() if i != 0]) for sample in df_samples]
            x_values = df_samples
            fig.add_trace(go.Bar(x=x_values, y=y_values, name=taxon, marker=dict(color=colors[c])))
            res.append(y_values)
    
    df_data = pd.DataFrame(res, index=all_taxa, columns=x_values)
    df_data_rel = df_data.div(df_data.sum(axis=0), axis=1) * 100

    # Specific Layout
    if fig_taxon == 'Hash':
        fig_taxon = st.session_state['TTT_hash'][:-1]
    fig.update_xaxes(dtick='linear')
    fig.update_yaxes(title=fig_mode)
    if fig_mode == 'Relative Reads':
        fig.update_layout(barnorm="percent")
    fig.update_layout(title='Distribution of Reads', barmode="stack")

    # Default Layout
    fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)

    # Plot
    st.plotly_chart(fig, config=st.session_state['TTT_config'])

    # Export
    export_plot("Read_proportions", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", fig, "plotly")
    export_table("Read_proportions", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", df_data_rel, "xlsx", True)

def read_hash_autocorrelation():

    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    TTT_colorscale = st.session_state['TTT_colorscale']
    TTT_colorsequence = st.session_state['TTT_colorsequence']
    TTT_color1 = st.session_state['TTT_color1']
    TTT_jitter = st.session_state['TTT_jitter']
    TTT_linewidth = st.session_state['TTT_linewidth']
    TTT_scattersize = st.session_state['TTT_scattersize'] * 0.5

    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table']
    df_samples = st.session_state['df_samples']
    TTT_hash = st.session_state['TTT_hash']

    # Per Sample
    x_values = []
    y_values = []
    for sample in df_samples:
        sub_df = df_taxon_table[sample]
        n_ESVs = len([i for i in sub_df.values.tolist() if i != 0])
        n_reads = sub_df.sum()
        x_values.append(n_reads)
        y_values.append(n_ESVs)

    fig = go.Figure()

    # Data
    fig.add_trace(go.Scatter(x=x_values, y=y_values, text=df_samples, mode='markers', cliponaxis=False, name='Samples', marker=dict(color=TTT_color1)))

    # Regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_values, y_values)
    x_line = np.array(x_values)
    y_line = intercept + slope * x_line
    fig.add_trace(go.Scatter(
        x=x_line,
        y=y_line,
        mode='lines',
        name='Regression line',
        line=dict(color=TTT_color2)))

    # Spearman
    rho, p = stats.spearmanr(x_values, y_values)
    rho, p = convert_r_p(rho, p)
    title = f'Read/{TTT_hash[:-1]} auto-correlation rho={rho} (p={p})'

    # Default Layout
    fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize), title=st.session_state['TTT_hash'])
    fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize), title='Reads')
    fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend, title=title)

    # Plot
    st.plotly_chart(fig, config=st.session_state['TTT_config'])

    # Export
    df = pd.DataFrame([df_samples, x_values, y_values]).transpose()
    df.columns = ['Sample', 'Reads', TTT_hash]
    export_plot("Read_proportions", f"{name}_reads_{TTT_hash}_corr", fig, "plotly")
    export_table("Read_proportions", f"{name}_reads_{TTT_hash}_corr", df, "xlsx", False)

def read_rarefaction():
    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    TTT_colorscale = st.session_state['TTT_colorscale']
    TTT_colorsequence = st.session_state['TTT_colorsequence']
    TTT_color1 = st.session_state['TTT_color1']
    TTT_jitter = st.session_state['TTT_jitter']
    TTT_linewidth = st.session_state['TTT_linewidth']
    TTT_scattersize = st.session_state['TTT_scattersize'] * 0.5

    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table']
    df_samples = st.session_state['df_samples']
    fig_taxon = st.session_state['read_rarefaction_taxon']
    read_rarefaction_splits = st.session_state['read_rarefaction_splits']
    read_rarefaction_reps = st.session_state['read_rarefaction_reps']
    fig_color = st.session_state['read_rarefaction_color']
    fig_display = st.session_state['read_rarefaction_display']

    if fig_taxon != 'Hash':
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)
    else:
        df_taxon_table_simple = df_taxon_table

    # Drop NaN
    df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']

    rarefaction_res = {}
    for sample in stqdm(df_samples, desc="Calculating Rarefaction Curves"):
        # Collect relevant data
        read_df = df_taxon_table_simple[[sample, fig_taxon]].copy()
        read_df = read_df[read_df[sample] != 0]
        max_taxa = read_df[fig_taxon].nunique()
        total_reads = read_df[sample].sum()

        # Expand reads into list for subsampling
        read_list = pd.Series(np.repeat(read_df[fig_taxon].to_list(), read_df[sample].to_list()))

        # Random Subsampling
        x_values = []
        y_values = []
        stdev_values = []
        sub_sample_size = total_reads / read_rarefaction_splits

        reads_pool = 0
        for subset in range(0, read_rarefaction_splits+1):
            x_res = []
            y_res = []
            for rep in range(0, read_rarefaction_reps):
                try:
                    sub_sample = read_list.sample(n=math.ceil(reads_pool))
                    n_taxa = len(set(sub_sample))
                except:
                    n_taxa = max_taxa
                if fig_display == 'Relative':
                    if max_taxa == 0:
                        n_taxa = 0
                    else:
                        n_taxa = n_taxa / max_taxa * 100
                x_res.append(math.ceil(reads_pool))
                y_res.append(n_taxa)
            # collect results
            avg_taxa = np.mean(y_res)
            stdev = np.std(y_res)
            x_values.append(math.ceil(reads_pool))
            y_values.append(avg_taxa)
            stdev_values.append(stdev)
            reads_pool += sub_sample_size

        df = pd.DataFrame([x_values, y_values, stdev_values]).transpose()
        df.columns = ['Reads', fig_taxon, 'STDEV']
        rarefaction_res[sample] = df

    # PLOT
    # Colors
    if fig_color == 'Color Scale':
        colors = get_colors_from_scale(TTT_colorscale, len(df_samples))
    if fig_color == 'Color Sequence':
        colors = get_colors_from_sequence(TTT_colorsequence, len(df_samples))
    if fig_color == 'Single Color':
        colors = [TTT_color1] * len(df_samples)

    fig = go.Figure()
    c=0
    for sample, df in rarefaction_res.items():
        x_values = df['Reads']
        y_values = df[fig_taxon]
        fig.add_trace(go.Scatter(x=x_values, y=y_values, name=sample, cliponaxis=False, mode='lines', line=dict(width=TTT_linewidth), marker=dict(color=colors[c])))
        c+=1

    if fig_display == 'Relative':
        fig.update_yaxes(range=(0,101), title=f'{fig_taxon} (%)')
    else:
        fig.update_yaxes(rangemode='tozero', title=fig_taxon)

    # Default Layout
    fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)

    # Plot
    st.plotly_chart(fig, config=st.session_state['TTT_config'])

    # Export
    export_plot("Read_proportions", f"{name}_{fig_taxon}_{fig_display.lower()}_read_rarefaction_curve", fig, "plotly")

    # Export multiple dfs to a single file
    try:
        active_project_path = st.session_state['active_project_path']
        export_file = active_project_path / "Read_proportions" / f"{name}_{fig_taxon}_{fig_display.lower()}_read_rarefaction_curve.xlsx"
        with pd.ExcelWriter(export_file, engine="xlsxwriter") as writer:
            i = 0
            for sample, df in rarefaction_res.items():
                sample_name = f"{i}_{sample[:28]}"
                df.to_excel(writer, sheet_name=str(sample_name), index=False)
                i += 1
    except:
        st.warning('Unable to write Rarefaction results to .xslx!')

########################################################################################################################
# Alpha Diversity

def alphadiv_diagram():

    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    TTT_colorscale = st.session_state['TTT_colorscale']
    TTT_colorsequence = st.session_state['TTT_colorsequence']
    TTT_color1 = st.session_state['TTT_color1']
    TTT_jitter = st.session_state['TTT_jitter']
    TTT_linewidth = st.session_state['TTT_linewidth']
    TTT_scattersize = st.session_state['TTT_scattersize'] * 0.5

    # Tables
    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    # Options
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples'].copy()
    fig_color = st.session_state['alphadiv_color'] # color sequence or scale
    fig_layout = st.session_state['alphadiv_layout'] # Bar or Box/Violin
    fig_level = st.session_state['alphadiv_level'] # sample or metadata
    fig_taxon = st.session_state['alphadiv_taxon'] # taxonomic level
    fig_nan = 'Exclude'
    fig_boxpoints = st.session_state['alphadiv_boxpoints']
    fig_measure = st.session_state['alphadiv_measure'] # richness or diff
    fig_mode = st.session_state['alphadiv_mode'] # quantification
    fig_interpration_lines = st.session_state['alphadiv_interpration']

    fig = go.Figure()

    # Handle Metadata
    df_taxon_table_simple = df_taxon_table.copy()
    if fig_nan == 'Exclude':
        df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']

    # calculate y-values
    all_taxa = sorted(df_taxon_table_simple[fig_taxon].unique())
    res = {}
    counts_dict = {}
    if fig_measure == 'Richness':
        for sample in df_samples:
            sub_df = df_taxon_table_simple[[fig_taxon, sample]]
            counts = sub_df[sub_df[sample] != 0][fig_taxon].nunique()
            res[sample] = counts
            counts = [len(sub_df[(sub_df[fig_taxon] == i) & (sub_df[sample] != 0)]) for i in all_taxa]
            counts_dict[sample] = counts # here we must provide an alternative df
    else:
        if fig_mode == 'Rel. Reads':
            for sample in df_samples:
                sub_df = df_taxon_table_simple[[fig_taxon, sample]]
                total_reads = sub_df[sample].sum()
                counts = [sub_df[sub_df[fig_taxon] == i][sample].sum() / total_reads * 100 for i in all_taxa]
                counts_dict[sample] = counts
                if fig_measure == 'Heip Evenness':
                    res[sample] = heip_e(counts)
                elif fig_measure == 'Pielou Evenness':
                    res[sample] = pielou_e(counts)
                else:
                    res[sample] = shannon(counts)
        else:
            # Haplotype Richness
            for sample in df_samples:
                sub_df = df_taxon_table_simple[[fig_taxon, sample]]
                counts = [len(sub_df[(sub_df[fig_taxon] == i) & (sub_df[sample] != 0)]) for i in all_taxa]
                counts_dict[sample] = counts
                if fig_measure == 'Heip Evenness':
                    res[sample] = heip_e(counts)
                elif fig_measure == 'Pielou Evenness':
                    res[sample] = pielou_e(counts)
                else:
                    res[sample] = shannon(counts)

    df1 = pd.DataFrame(counts_dict.values(), index=counts_dict.keys(), columns=all_taxa).transpose()
    df2 = pd.DataFrame(res.values(), index=res.keys(), columns=[fig_measure])

    if fig_layout == 'Bar':
        # Colors
        x_values = list(res.keys())
        if fig_color == 'Color Scale':
            colors = get_colors_from_scale(TTT_colorscale, len(x_values))
        if fig_color == 'Color Sequence':
            colors = get_colors_from_sequence(TTT_colorsequence, len(x_values))
        if fig_color == 'Single Color':
            colors = [TTT_color1] * len(x_values)

        # Traces
        c = 0
        for x_values, y_values in res.items():
            fig.add_trace(go.Bar(x=[x_values], y=[y_values], name=x_values, marker=dict(color=colors[c])))
            c += 1

        # Specific Layout
        if fig_taxon == 'Hash':
            fig_taxon = st.session_state['TTT_hash'][:-1]
            fig.update_yaxes(title=f"{fig_taxon} {fig_measure}")
        else:
            fig.update_yaxes(title=f"{st.session_state['tax_level_plurals'][fig_taxon]} {fig_measure}")
        if fig_measure == 'Heip Evenness' or fig_measure == 'Pielou Evenness':
            fig.update_yaxes(range=(0,1.1))
        else:
            fig.update_yaxes(rangemode='tozero')
        fig.update_xaxes(dtick='linear')

        # Add intrepration guidance
        if fig_interpration_lines == "Yes" and (fig_measure == 'Heip Evenness' or fig_measure == 'Pielou Evenness'):
            fig.add_hrect(y0=0.8, y1=1.0, fillcolor="green", opacity=0.2, layer="below", line_width=0)
            fig.add_hrect(y0=0.5, y1=0.8, fillcolor="yellow", opacity=0.2, layer="below", line_width=0)
            fig.add_hrect(y0=0.0, y1=0.5, fillcolor="orange", opacity=0.2, layer="below", line_width=0)

        # Default Layout
        fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)

        # Plot
        st.plotly_chart(fig, config=st.session_state['TTT_config'])

        # Export
        export_plot("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}", fig, "plotly")
        export_table("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}_1", df1, "xlsx", True)
        export_table("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}_2", df2, "xlsx", True)

    elif fig_layout == 'Scatter':
        # Colors
        x_values = list(res.keys())
        y_values = pd.Series(res).fillna(0).tolist()
        fig.add_trace(go.Scatter(x=x_values, y=y_values, fill='tozeroy', name='', line=dict(width=TTT_linewidth), marker=dict(color=TTT_color1)))

        # Specific Layout
        if fig_taxon == 'Hash':
            fig_taxon = st.session_state['TTT_hash'][:-1]
            fig.update_yaxes(title=f"{fig_taxon} {fig_measure}")
        else:
            fig.update_yaxes(title=f"{st.session_state['tax_level_plurals'][fig_taxon]} {fig_measure}")
        if fig_measure == 'Heip Evenness':
            fig.update_yaxes(range=(0,1.1))
        else:
            fig.update_yaxes(rangemode='tozero')
        fig.update_xaxes(dtick='linear')

        # Add intrepration guidance
        if fig_interpration_lines == "Yes" and (fig_measure == 'Heip Evenness' or fig_measure == 'Pielou Evenness'):
            fig.add_hrect(y0=0.8, y1=1.0, fillcolor="green", opacity=0.2, layer="below", line_width=0)
            fig.add_hrect(y0=0.5, y1=0.8, fillcolor="yellow", opacity=0.2, layer="below", line_width=0)
            fig.add_hrect(y0=0.0, y1=0.5, fillcolor="orange", opacity=0.2, layer="below", line_width=0)

        # Default Layout
        fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_layout(template=template, font_size=fontsize, showlegend=False)

        # Plot
        st.plotly_chart(fig, config=st.session_state['TTT_config'])

        # Export
        export_plot("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}", fig, "plotly")
        export_table("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}_1", df1, "xlsx", True)
        export_table("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}_2", df2, "xlsx", True)

    elif fig_layout == 'Box' or fig_layout == 'Violin':
        if fig_level == 'Samples':
            st.error('Please select a metadata category instead of "Samples"!')
            return

        # Collect metadata
        test_groups = [i for i in list(df_metadata_table[fig_level].unique()) if i != '']
        # Colors
        if fig_color == 'Color Scale':
            colors = get_colors_from_scale(TTT_colorscale, len(test_groups))
        if fig_color == 'Color Sequence':
            colors = get_colors_from_sequence(TTT_colorsequence, len(test_groups))
        if fig_color == 'Single Color':
            colors = [TTT_color1] * len(test_groups)

        c = 0
        for test_group in test_groups:
            group_samples = df_metadata_table[df_metadata_table[fig_level] == test_group]['Samples'].values.tolist()
            y_values = df2[df2.index.isin(group_samples)][fig_measure].values.tolist()
            if fig_layout == 'Box':
                boxpoints = 'all' if fig_boxpoints == 'All' else False
                fig.add_trace(go.Box(y=y_values, name=test_group, jitter=TTT_jitter, boxpoints=boxpoints, line=dict(width=TTT_linewidth), marker=dict(size=TTT_scattersize, color=colors[c])))
            else:
                boxpoints = 'all' if fig_boxpoints == 'All' else False
                fig.add_trace(go.Violin(y=y_values, name=test_group, jitter=TTT_jitter, points=boxpoints, line=dict(width=TTT_linewidth), marker=dict(size=TTT_scattersize, color=colors[c])))
            c += 1

        # Specific Layout
        if fig_taxon == 'Hash':
            fig_taxon = st.session_state['TTT_hash'][:-1]
            fig.update_yaxes(title=f"{fig_taxon} {fig_measure}")
        else:
            fig.update_yaxes(title=f"{st.session_state['tax_level_plurals'][fig_taxon]} {fig_measure}")
        if fig_measure == 'Heip Evenness':
            fig.update_yaxes(range=(0,1.05))
        else:
            fig.update_yaxes(rangemode='tozero')
        fig.update_xaxes(dtick='linear')

        # Add intrepration guidance
        if fig_interpration_lines == "Yes" and (fig_measure == 'Heip Evenness' or fig_measure == 'Pielou Evenness'):
            fig.add_hrect(y0=0.8, y1=1.0, fillcolor="green", opacity=0.2, layer="below", line_width=0)
            fig.add_hrect(y0=0.5, y1=0.8, fillcolor="yellow", opacity=0.2, layer="below", line_width=0)
            fig.add_hrect(y0=0.0, y1=0.5, fillcolor="orange", opacity=0.2, layer="below", line_width=0)

        # Default Layout
        fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)

        # Plot
        st.plotly_chart(fig, config=st.session_state['TTT_config'])

        # Export
        export_plot("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}", fig, "plotly")
        export_table("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}_1", df1, "xlsx", True)
        export_table("Alpha_diversity", f"{name}_{fig_taxon}_{fig_measure}_{fig_mode}_2", df2, "xlsx", True)

def pycircos_plot():

    # Layout
    fontsize = st.session_state['TTT_fontsize']

    # Tables
    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_samples = st.session_state['df_samples']


    # Options
    fig_taxon0 = st.session_state['pycircos_taxon0'] # Outer Groups
    fig_taxon1 = st.session_state['pycircos_taxon1'] # Inner Groups
    fig_taxon2 = st.session_state['pycircos_taxon2'] # Diversity Measure
    fig_nan = st.session_state['pycircos_nan'] # include or exclude
    fig_colors = st.session_state['pycircos_colors']

    if fig_taxon2 != 'Hash':
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)
    else:
        df_taxon_table_simple = df_taxon_table.copy()

    # Analysis
    if fig_nan == 'Exclude':
        df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon2] != '']

    # Test groups
    test_groups = [i for i in df_taxon_table_simple[fig_taxon0].unique() if i != '']

    res = {}
    ## loop through all groups
    for group in test_groups:
        # Collect data
        group_res = []
        sub_df = df_taxon_table_simple[df_taxon_table_simple[fig_taxon0] == group]
        taxa = [i for i in sub_df[fig_taxon1].unique() if i != '']
        if taxa == []:
            continue
        total_reads = sub_df[df_samples].sum().sum()

        for taxon in taxa:
            taxon_df = sub_df[sub_df[fig_taxon1] == taxon]
            richness = len([i for i in taxon_df[fig_taxon2].unique() if i != ''])
            n_reads = taxon_df[df_samples].sum().sum()
            rel_reads = n_reads / total_reads * 100
            group_res.append([taxon, richness, rel_reads])
        res_df = pd.DataFrame(group_res, columns=['Taxon', 'Richness', 'Reads'])
        res[group] = res_df

    sectors = {key:len(values) for key,values in res.items()}
    circos = Circos(sectors, space=1)
    for sector, group in zip(circos.sectors, res.keys()):
        # color
        if fig_colors == 'custom':
            color = st.session_state[f'pycircos_colors_{group}']
            cmap = st.session_state['pycircos_cmaps_dict'][color]
        else:
            color = st.session_state['pycircos_colors']
            cmap = st.session_state['pycircos_cmaps_dict'][color]

        # arrange sector
        x = np.arange(sector.start, sector.end) + 0.5

        # Add Species heatmap
        y = res[sector.name]['Richness'].values.tolist()
        t = res[sector.name]['Taxon'].values.tolist()
        matrix = np.array([y])
        r = 51
        heatmap_track = sector.add_track((r, r+20))
        heatmap_track.axis()
        vmax = max(max(matrix))
        heatmap_track.heatmap(matrix, cmap=cmap, vmin=0, vmax=vmax)
        heatmap_track.xticks(x, labels=t, label_orientation='vertical', label_size=fontsize)

        # Add labels in between  -> Species
        y = res[sector.name]['Richness'].values.tolist()
        heatmap_track = sector.add_track((r, r))
        heatmap_track.axis()
        heatmap_track.bar(x, y, color=color, vmin=0, vmax=max(y))
        heatmap_track.xticks(x, labels=[str(i) for i in y], label_orientation='vertical', label_size=fontsize-2)

        # Add Reads heatmap
        y = res[sector.name]['Reads'].values.tolist()
        matrix = np.array([y])
        r = 28
        heatmap_track = sector.add_track((r, r+20))
        heatmap_track.axis()
        heatmap_track.heatmap(matrix, cmap=cmap, vmin=0, vmax=100)

        # Add labels in between -> Reads
        y = [round(i,2) for i in res[sector.name]['Reads'].values.tolist()]
        labels = [str(i) if i >= 0.01 else "<0.01" for i in y]
        heatmap_track = sector.add_track((r, r))
        heatmap_track.axis()
        heatmap_track.bar(x, y, color=color, vmax=100)
        heatmap_track.xticks(x, labels=labels, label_orientation='vertical', label_size=fontsize-4)

        # Plot inside track
        r = 21
        track = sector.add_track((r, r+5))
        if st.session_state['pycircos_inner_labels'] == 'Show':
            track.text(sector.name, size=9)
        track.axis(fc=color, ec="black", alpha=0.5)

    fig = circos.plotfig(figsize=(8, 8))
    st.pyplot(fig)

    # Export
    export_plot("Alpha_diversity", f"{name}_circos_{fig_taxon0}_{fig_taxon1}_{fig_taxon2}", fig, "matplot")

########################################################################################################################
# Beta Diversity

def beta_diversity_matrix():

    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']

    # Tables
    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    # Options
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples'].copy()
    fig_labels = st.session_state['betadiv_labels'] # color sequence or scale
    fig_mode = st.session_state['betadiv_mode'] # jaccard or bray curtis
    fig_measure = st.session_state['betadiv_measure'] # spearman or bray curtis
    fig_level = st.session_state['betadiv_level'] # sample or metadata
    fig_taxon = st.session_state['betadiv_taxon'] # taxonomic level
    fig_nan = st.session_state['betadiv_nan'] # include or exclude

    # Handle Metadata
    if fig_level != 'Samples':
        df_taxon_table, df_samples = concatenate_by_metadata(fig_level)

    if fig_taxon == 'Hash':
        df_taxon_table_simple = df_taxon_table.copy()
    else:
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)

    # Analysis
    if fig_nan == 'Exclude':
        df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']

    df_taxon_table_simple_s = (df_taxon_table_simple[[fig_taxon] + df_samples].groupby(fig_taxon, as_index=False).sum())
    richness_list = {i:df_taxon_table_simple_s[i].values.tolist() for i in df_samples}

    matrix = []
    if fig_mode == 'Jaccard (qualitative)':
        for sample1 in df_samples:
            row = []
            v1 = [1 if i !=0 else 0 for i in richness_list[sample1]]
            for sample2 in df_samples:
                v2 = [1 if i !=0 else 0 for i in richness_list[sample2]]
                if fig_measure == 'Similarity':
                    d = 1-distance.jaccard(v1, v2)
                else:
                    d = distance.jaccard(v1, v2)
                row.append(d)
            matrix.append(row)
    else:
        for sample1 in df_samples:
            row = []
            v1 = richness_list[sample1]
            for sample2 in df_samples:
                v2 = richness_list[sample2]
                if fig_measure == 'Similarity':
                    d = 1-distance.jaccard(v1, v2)
                else:
                    d = distance.braycurtis(v1, v2)
                row.append(d)
            matrix.append(row)

    # Output Table
    df = pd.DataFrame(matrix, columns=df_samples, index=df_samples)
    df = df.round(2)

    # Output Plot
    text_auto = True if fig_labels == 'Show' else False
    fig = px.imshow(df, text_auto=text_auto, color_continuous_scale=st.session_state['TTT_colorscale'])

    # Specific Layout
    if fig_taxon == 'Hash':
        fig_taxon = st.session_state['TTT_hash'][:-1]
    fig.update_xaxes(dtick='linear')
    fig.update_yaxes(dtick='linear')

    # Default Layout
    fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)

    # Plot
    st.plotly_chart(fig, config=st.session_state['TTT_config'])

    # Export
    export_plot("Beta_diversity", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", fig, "plotly")
    export_table("Beta_diversity", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", df, "xlsx", True)

def draw_outlines(fig, x_values, y_values, metadata, color):
    try:
        ## collect samples that form the outline
        x_plane, y_plane = [], []
        hull = ConvexHull(np.column_stack((x_values, y_values)))
        for i in hull.vertices:
            x_plane.append(x_values[i])
            y_plane.append(y_values[i])
        x_plane.append(x_values[hull.vertices[0]])
        y_plane.append(y_values[hull.vertices[0]])

        ## draw the outline
        fig.add_trace(go.Scatter(x=x_plane, y=y_plane, mode='lines', name=metadata, marker_color=color, fill='toself'))
    except:
        pass

def pcoa_calculation():

    # Tables
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    # Options
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples'].copy()
    fig_mode = st.session_state['pcoa_mode'] # spearman or bray curtis
    fig_taxon = st.session_state['pcoa_taxon'] # taxonomic level
    fig_level = st.session_state['pcoa_level']

    # Handle Metadata
    if fig_level != 'Samples':
        drop_samples = df_metadata_table[df_metadata_table[fig_level] == '']['Samples'].values.tolist()
        if len(drop_samples) != 0:
            df_samples = [i for i in df_samples if i not in drop_samples]
            df_taxon_table = df_taxon_table.drop(columns=drop_samples, errors='ignore')

    if fig_taxon == 'Hash':
        df_taxon_table_simple = df_taxon_table.copy()
    else:
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)
    df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']
    df_taxon_table_simple_s = (df_taxon_table_simple[[fig_taxon] + df_samples].groupby(fig_taxon, as_index=False).sum())
    richness_list = {i:df_taxon_table_simple_s[i].values.tolist() for i in df_samples}

    matrix = []
    if fig_mode == 'Jaccard (qualitative)':
        for sample1 in df_samples:
            row = []
            v1 = [1 if i !=0 else 0 for i in richness_list[sample1]]
            for sample2 in df_samples:
                v2 = [1 if i !=0 else 0 for i in richness_list[sample2]]
                d = distance.jaccard(v1, v2)
                row.append(d)
            matrix.append(row)
    else:
        for sample1 in df_samples:
            row = []
            v1 = richness_list[sample1]
            for sample2 in df_samples:
                v2 = richness_list[sample2]
                d = distance.braycurtis(v1, v2)
                row.append(d)
            matrix.append(row)

    # Output Table
    distance_matrix_df = pd.DataFrame(matrix, columns=df_samples, index=df_samples)

    # Calculate pcoa
    pcoa_res = pcoa(distance_matrix_df)

    # Collect expained variance
    pcoa_explained_variance_df = pd.DataFrame(pcoa_res.proportion_explained, columns=['explained_variance']) * 100
    pcoa_explained_variance_df = pcoa_explained_variance_df[pcoa_explained_variance_df['explained_variance'] > 1]

    # Collect values and already filter for PC axes with >1% explained variance
    pcoa_df = pd.DataFrame(pcoa_res.samples)
    pcoa_df.index = df_samples
    pcoa_axes = list(pcoa_explained_variance_df.index)
    pcoa_df = pcoa_df[pcoa_axes]
    explained_variance_df = pd.DataFrame(pcoa_explained_variance_df)

    return pcoa_df, explained_variance_df, distance_matrix_df

def pcoa_plot():

    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    scattersize = st.session_state['TTT_scattersize']
    showlabels = st.session_state['pcoa_showlabels']
    linewidth = st.session_state['TTT_linewidth']

    # Extract axes to display
    name = st.session_state['table_display_name'].stem
    fig_mode = st.session_state['pcoa_mode'] # spearman or bray curtis
    fig_taxon = st.session_state['pcoa_taxon'] # taxonomic level
    pcoa_df = st.session_state['pcoa_df']
    df_samples = pcoa_df.index.tolist()
    pcoa_explained_variance_dict = st.session_state['pcoa_explained_variance_dict']
    pcoa_distance_matrix_df = st.session_state['pcoa_distance_matrix_df']
    pcoa_x = st.session_state['pcoa_x']
    pcoa_y = st.session_state['pcoa_y']
    pcoa_z = st.session_state['pcoa_z']
    pcoa_level = st.session_state['pcoa_level']
    df_metadata_table = st.session_state['df_metadata_table']
    sample_categories = [i for i in df_metadata_table[pcoa_level].values.tolist() if i != '']

    # ANOSIM
    if pcoa_level != 'Samples':
        dm = DistanceMatrix(pcoa_distance_matrix_df)
        anosim_result = anosim(dm, sample_categories, permutations=999)
    else:
        anosim_result = ''

    ## Display anosim
    try:
        r_value = round(anosim_result['test statistic'], 3)
        p_value = anosim_result['p-value']
        title = f'{pcoa_level}: ANOSIM R={r_value}, p={p_value}'
    except:
        title = f'ANOSIM not calculated'

    color_mode = st.session_state['pcoa_color']
    if pcoa_level == 'Samples':
        color_mode = 'Single Color'
    if color_mode == 'Color scale':
        colors = get_colors_from_scale(st.session_state['TTT_colorscale'], len(set(sample_categories)))
        pcoa_df['Group'] = [i for i in df_metadata_table[pcoa_level].values.tolist() if i != '']
        group_colors = {j: colors[i] for i, j in enumerate(set(sample_categories))}
        pcoa_df['Colors'] = [group_colors[i] for i in pcoa_df['Group'].values]
    elif color_mode == 'Color sequence':
        colors = get_colors_from_sequence(st.session_state['TTT_colorsequence'], len(set(sample_categories)))
        pcoa_df['Group'] = [i for i in df_metadata_table[pcoa_level].values.tolist() if i != '']
        group_colors = {j: colors[i] for i, j in enumerate(set(sample_categories))}
        pcoa_df['Colors'] = [group_colors[i] for i in pcoa_df['Group'].values]
    else:
        colors = [st.session_state['TTT_color1']] * len(df_samples)
        pcoa_df['Colors'] = colors
        pcoa_df['Group'] = [i for i in df_metadata_table[pcoa_level].values.tolist() if i != '']

    if pcoa_z == 'None':
        x_title = pcoa_explained_variance_dict[pcoa_x]
        y_title = pcoa_explained_variance_dict[pcoa_y]
        pcoa_2D_df = pcoa_df[[x_title, y_title, 'Colors', 'Group']]

        fig = go.Figure()

        for metadata in pcoa_df['Group'].unique():
            sub_df = pcoa_2D_df.loc[pcoa_2D_df['Group'] == metadata]
            x_values = sub_df[x_title].values.tolist()
            y_values = sub_df[y_title].values.tolist()
            text_values = list(sub_df.index)
            colors = sub_df['Colors'].values.tolist()

            fig.add_trace(go.Scatter(x=x_values, y=y_values, text=text_values, marker_color=colors, marker=dict(size=scattersize), name=metadata, mode='markers'))

        if pcoa_level != 'Samples' and pcoa_groups == 'Outlines':
            for category in set(pcoa_df['Group'].values.tolist()):
                sub_df = pcoa_df.loc[pcoa_df['Group'] == category]
                if len(sub_df) >= 2:
                    x_values = sub_df[x_title].values.tolist()
                    y_values = sub_df[y_title].values.tolist()
                    color = sub_df['Colors'].drop_duplicates().values.tolist()[0]
                    draw_outlines(fig, x_values, y_values, category, color)

        elif pcoa_level != 'Samples' and pcoa_groups == 'Continous':
            x_values = pcoa_df[x_title].values.tolist()
            y_values = pcoa_df[y_title].values.tolist()
            fig.add_trace(go.Scatter(x=x_values, y=y_values, marker=dict(color=st.session_state['TTT_color1']), line=dict(width=linewidth), mode='lines'))

        if showlabels != 'Hide':
            x_values = pcoa_df[x_title].values.tolist()
            y_values = pcoa_df[y_title].values.tolist()
            t_values = list(pcoa_df.index)
            fig.add_trace(go.Scatter(x=x_values, y=y_values, text=t_values, textposition=showlabels, mode='text'))

        # Update layout
        fig.update_xaxes(title=pcoa_x)
        fig.update_yaxes(title=pcoa_y)
        fig.update_layout(title=title,
                          template=template,
                          font_size=fontsize,
                          showlegend=showlegend,
                          yaxis_title_font=dict(size=fontsize),
                          yaxis_tickfont=dict(size=fontsize),
                          xaxis_title_font=dict(size=fontsize),
                          xaxis_tickfont=dict(size=fontsize))

        # Plot
        st.plotly_chart(fig, config=st.session_state['TTT_config'])
        export_plot("Beta_diversity", f"{name}_{fig_mode}_{fig_taxon}", fig, "plotly")
        export_table("Beta_diversity", f"{name}_{fig_mode}_{fig_taxon}", pcoa_df, "xlsx")

def nmds_calculation():

    # Tables
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    # Options
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples'].copy()
    fig_mode = st.session_state['nmds_mode'] # spearman or bray curtis
    fig_taxon = st.session_state['nmds_taxon'] # taxonomic level
    fig_level = st.session_state['nmds_level']
    nmds_axes = st.session_state['nmds_axes']

    # Handle Metadata
    if fig_level != 'Samples':
        drop_samples = df_metadata_table[df_metadata_table[fig_level] == '']['Samples'].values.tolist()
        if len(drop_samples) != 0:
            df_samples = [i for i in df_samples if i not in drop_samples]
            df_taxon_table = df_taxon_table.drop(columns=drop_samples, errors='ignore')

    if fig_taxon == 'Hash':
        df_taxon_table_simple = df_taxon_table.copy()
    else:
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)
    df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']
    df_taxon_table_simple_s = (df_taxon_table_simple[[fig_taxon] + df_samples].groupby(fig_taxon, as_index=False).sum())
    richness_list = {i:df_taxon_table_simple_s[i].values.tolist() for i in df_samples}

    matrix = []
    if fig_mode == 'Jaccard (qualitative)':
        for sample1 in df_samples:
            row = []
            v1 = [1 if i !=0 else 0 for i in richness_list[sample1]]
            for sample2 in df_samples:
                v2 = [1 if i !=0 else 0 for i in richness_list[sample2]]
                d = distance.jaccard(v1, v2)
                row.append(d)
            matrix.append(row)
    else:
        for sample1 in df_samples:
            row = []
            v1 = richness_list[sample1]
            for sample2 in df_samples:
                v2 = richness_list[sample2]
                d = distance.braycurtis(v1, v2)
                row.append(d)
            matrix.append(row)

    # Output Table
    distance_matrix_df = pd.DataFrame(matrix, columns=df_samples, index=df_samples)

    # Calculate pcoa
    nmds = MDS(
        n_components=nmds_axes,
        metric=False,
        dissimilarity="precomputed",
        random_state=42,
        n_init=10,
        max_iter=3000)

    coords = nmds.fit_transform(distance_matrix_df.values)
    if nmds_axes == 2:
        nmds_df = pd.DataFrame(coords, index=df_samples, columns=["NMDS1", "NMDS2"])
    else:
        nmds_df = pd.DataFrame(coords, index=df_samples, columns=["NMDS1", "NMDS2", "NMDS3"])
    nmds_stress = nmds.stress_

    return nmds_df, nmds_stress

def nmds_plot():

    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    scattersize = st.session_state['TTT_scattersize']
    showlabels = st.session_state['nmds_showlabels']
    linewidth = st.session_state['TTT_linewidth']
    nmds_axes = st.session_state['nmds_axes']

    # Extract axes to display
    name = st.session_state['table_display_name'].stem
    fig_mode = st.session_state['nmds_mode'] # spearman or bray curtis
    fig_taxon = st.session_state['nmds_taxon'] # taxonomic level
    nmds_level = st.session_state['nmds_level']
    nmds_df = st.session_state['nmds_df']
    nmds_stress = st.session_state['nmds_stress']

    df_samples = nmds_df.index.tolist()
    df_metadata_table = st.session_state['df_metadata_table']
    sample_categories = [i for i in df_metadata_table[nmds_level].values.tolist() if i != '']

    # Create a title
    title = f'Stress: {nmds_stress:.3}'

    color_mode = st.session_state['nmds_color']
    if nmds_level == 'Samples':
        color_mode = 'Single Color'
    if color_mode == 'Color scale':
        colors = get_colors_from_scale(st.session_state['TTT_colorscale'], len(set(sample_categories)))
        nmds_df['Group'] = [i for i in df_metadata_table[nmds_level].values.tolist() if i != '']
        group_colors = {j: colors[i] for i, j in enumerate(set(sample_categories))}
        nmds_df['Colors'] = [group_colors[i] for i in nmds_df['Group'].values]
    elif color_mode == 'Color sequence':
        colors = get_colors_from_sequence(st.session_state['TTT_colorsequence'], len(set(sample_categories)))
        nmds_df['Group'] = [i for i in df_metadata_table[nmds_level].values.tolist() if i != '']
        group_colors = {j: colors[i] for i, j in enumerate(set(sample_categories))}
        nmds_df['Colors'] = [group_colors[i] for i in nmds_df['Group'].values]
    else:
        colors = [st.session_state['TTT_color1']] * len(df_samples)
        nmds_df['Colors'] = colors
        nmds_df['Group'] = [i for i in df_metadata_table[nmds_level].values.tolist() if i != '']

    if nmds_axes == 2:
        x_title = 'NMDS1'
        y_title = 'NMDS2'
        nmds_2D_df = nmds_df[[x_title, y_title, 'Colors', 'Group']]

        fig = go.Figure()

        for metadata in nmds_df['Group'].unique():
            sub_df = nmds_2D_df.loc[nmds_2D_df['Group'] == metadata]
            x_values = sub_df[x_title].values.tolist()
            y_values = sub_df[y_title].values.tolist()
            text_values = list(sub_df.index)
            colors = sub_df['Colors'].values.tolist()

            fig.add_trace(go.Scatter(x=x_values, y=y_values, text=text_values, marker_color=colors, marker=dict(size=scattersize), name=metadata, mode='markers'))

        if nmds_level != 'Samples' and nmds_groups == 'Outlines':
            for category in set(nmds_df['Group'].values.tolist()):
                sub_df = nmds_df.loc[nmds_df['Group'] == category]
                if len(sub_df) >= 2:
                    x_values = sub_df[x_title].values.tolist()
                    y_values = sub_df[y_title].values.tolist()
                    color = sub_df['Colors'].drop_duplicates().values.tolist()[0]
                    draw_outlines(fig, x_values, y_values, category, color)

        elif nmds_level != 'Samples' and nmds_groups == 'Continous':
            x_values = nmds_df[x_title].values.tolist()
            y_values = nmds_df[y_title].values.tolist()
            fig.add_trace(go.Scatter(x=x_values, y=y_values, marker=dict(color=st.session_state['TTT_color1']), line=dict(width=linewidth), mode='lines'))

        if showlabels != 'Hide':
            x_values = nmds_df[x_title].values.tolist()
            y_values = nmds_df[y_title].values.tolist()
            t_values = list(nmds_df.index)
            fig.add_trace(go.Scatter(x=x_values, y=y_values, text=t_values, textposition=showlabels, mode='text'))

        # Update layout
        fig.update_xaxes(title=x_title)
        fig.update_yaxes(title=y_title)
        fig.update_layout(title=title,
                          template=template,
                          font_size=fontsize,
                          showlegend=showlegend,
                          yaxis_title_font=dict(size=fontsize),
                          yaxis_tickfont=dict(size=fontsize),
                          xaxis_title_font=dict(size=fontsize),
                          xaxis_tickfont=dict(size=fontsize))

        # Plot
        st.plotly_chart(fig, config=st.session_state['TTT_config'])
        export_plot("Beta_diversity", f"{name}_{fig_mode}_{fig_taxon}", fig, "plotly")
        export_table("Beta_diversity", f"{name}_{fig_mode}_{fig_taxon}", nmds_df, "xlsx")

    else:
        x_title = 'NMDS1'
        y_title = 'NMDS2'
        z_title = 'NMDS3'
        nmds_2D_df = nmds_df[[x_title, y_title, z_title, 'Colors', 'Group']]

        fig = go.Figure()

        for metadata in nmds_df['Group'].unique():
            sub_df = nmds_2D_df.loc[nmds_2D_df['Group'] == metadata]
            x_values = sub_df[x_title].values.tolist()
            y_values = sub_df[y_title].values.tolist()
            z_values = sub_df[z_title].values.tolist()
            text_values = list(sub_df.index)
            colors = sub_df['Colors'].values.tolist()

            fig.add_trace(go.Scatter3d(x=x_values, y=y_values, z=z_values, text=text_values, marker_color=colors, marker=dict(size=scattersize), name=metadata, mode='markers'))

        # Update layout
        fig.update_xaxes(title=x_title)
        fig.update_yaxes(title=y_title)
        fig.update_layout(title=title,
                          template=template,
                          font_size=fontsize,
                          showlegend=showlegend,
                          yaxis_title_font=dict(size=fontsize),
                          yaxis_tickfont=dict(size=fontsize),
                          xaxis_title_font=dict(size=fontsize),
                          xaxis_tickfont=dict(size=fontsize))

        # Plot
        st.plotly_chart(fig, config=st.session_state['TTT_config'])
        export_plot("Beta_diversity", f"{name}_{fig_mode}_{fig_taxon}", fig, "plotly")
        export_table("Beta_diversity", f"{name}_{fig_mode}_{fig_taxon}", nmds_df, "xlsx")

########################################################################################################################
# Population Dynamics
def haplotype_distribution():
    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    TTT_colorscale = st.session_state['TTT_colorscale']
    TTT_colorsequence = st.session_state['TTT_colorsequence']
    TTT_color1 = st.session_state['TTT_color1']
    TTT_jitter = st.session_state['TTT_jitter']
    TTT_linewidth = st.session_state['TTT_linewidth']
    TTT_scattersize = st.session_state['TTT_scattersize'] * 0.5
    TTT_hash = st.session_state['TTT_hash']

    # Tables
    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    # Options
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples'].copy()
    haplotype_taxa = st.session_state['haplotype_taxa']
    haplotype_level = st.session_state['haplotype_level']
    name = st.session_state['table_display_name'].stem

    colors = get_colors_from_scale(TTT_colorscale, len(haplotype_taxa))

    fig = go.Figure()
    output_lst = []
    for c, haplotype_taxon in enumerate(haplotype_taxa):
        sub_df = df_taxon_table[df_taxon_table[haplotype_level] == haplotype_taxon]

        y_values = []
        for sample in df_samples:
            y_value = len([i for i in sub_df[sample].values.tolist() if i != 0])
            y_values.append(y_value)


        if haplotype_mode == 'Scatter':
            fig.add_trace(go.Scatter(x=df_samples, y=y_values,
                                     marker=dict(color=colors[c]),
                                     name=haplotype_taxon, cliponaxis=False))
            fig.update_yaxes(rangemode='tozero', title=TTT_hash)
            fig.update_xaxes(dtick='linear')

        elif haplotype_mode == 'Bar':
            fig.add_trace(go.Bar(x=df_samples, y=y_values,
                                 marker=dict(color=colors[c]),
                                 name=haplotype_taxon))
            fig.update_yaxes(rangemode='tozero', title=TTT_hash)
            fig.update_xaxes(dtick='linear')
            fig.update_layout(barmode='stack')
        output_lst.append([haplotype_taxon] + y_values)

    # Default Layout
    fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
    fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend, title=f'{haplotype_level}-level {TTT_hash}')

    # Plot
    st.plotly_chart(fig, config=st.session_state['TTT_config'])
    export_plot('Haplotypes', f'{name}_{haplotype_level}_{haplotype_mode}', fig, 'plotly')

    # Dataframe
    res_df = pd.DataFrame(output_lst, columns=['Taxon'] + df_samples)
    export_table("Haplotypes", f'{name}_{haplotype_level}_{haplotype_mode}', res_df, "xlsx")


########################################################################################################################
# Time Series
def time_series_richness():

    # Layout
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    scattersize = st.session_state['TTT_scattersize']
    linewidth = st.session_state['TTT_linewidth']
    scattersize = st.session_state['TTT_scattersize']
    color1 = st.session_state['TTT_color1']
    color2 = st.session_state['TTT_color2']

    # Tables
    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_samples = st.session_state['df_samples']

    # Options
    fig_level = st.session_state['ts_level']
    fig_taxon = st.session_state['ts_taxon']
    fig_confidence_interval = st.session_state['ts_ci']
    fig_bootstraps = st.session_state['ts_bootstraps']

    # Handle Metadata
    if fig_level != 'Samples':
        df_taxon_table, df_samples = concatenate_by_metadata(fig_level)

    if fig_taxon != 'Hash':
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)
    else:
        df_taxon_table_simple = df_taxon_table.copy()

    # Remove all 'nan'
    df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']

    # Prepare x-values and y-values (using actual sample names as x-values)
    x_values = []
    y_values = []
    for sample in df_samples:
        sub_df = df_taxon_table_simple[[fig_taxon, sample]]
        measure = len(sub_df[sub_df[sample] != 0].values.tolist())
        x_values.append(sample)
        y_values.append(measure)

    # Plot original data points
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_values, y=y_values, name='Original Data', mode='markers', marker=dict(color=color1, size=scattersize)))

    # Calculate LOESS (Locally Estimated Scatterplot Smoothing) for the original data
    # Since x_values are now sample names, we need numerical indices for LOESS smoothing
    x_numeric = list(range(len(x_values)))
    loess_smoothed = sm.nonparametric.lowess(endog=y_values, exog=x_numeric, frac=0.3)
    loess_x_numeric = loess_smoothed[:, 0]
    loess_y = loess_smoothed[:, 1]

    # Map numerical x-values back to sample names
    loess_x = [x_values[int(i)] for i in loess_x_numeric]

    # Add LOESS-smoothed line to the plot
    fig.add_trace(go.Scatter(x=loess_x, y=loess_y, name='LOESS Smoothing', mode='lines', line=dict(color=color2, width=linewidth)))

    # Bootstrapping to calculate confidence intervals
    bootstrapped_fits = []
    for _ in range(fig_confidence_interval):
        # Resample the data with replacement
        resample_indices = np.random.choice(range(len(x_numeric)), size=len(x_numeric), replace=True)
        resample_x_numeric = np.array(x_numeric)[resample_indices]
        resample_y = np.array(y_values)[resample_indices]

        # Fit LOESS on resampled data
        loess_bootstrap = sm.nonparametric.lowess(endog=resample_y, exog=resample_x_numeric, frac=0.3)

        # Interpolate to original x-numeric values
        interp_bootstrap_y = np.interp(x_numeric, loess_bootstrap[:, 0], loess_bootstrap[:, 1])
        bootstrapped_fits.append(interp_bootstrap_y)

    # Convert list to numpy array for easier calculations
    bootstrapped_fits = np.array(bootstrapped_fits)

    # Calculate the confidence intervals (CI) at each x-value
    lower_bound = np.percentile(bootstrapped_fits, (100 - fig_confidence_interval) / 2, axis=0)
    upper_bound = np.percentile(bootstrapped_fits, 100 - (100 - fig_confidence_interval) / 2, axis=0)

    # Add shaded region for the confidence interval
    fig.add_trace(go.Scatter(
        x=x_values + x_values[::-1],  # x values forward and backward for the fill area
        y=np.concatenate([upper_bound, lower_bound[::-1]]),  # upper bound followed by reversed lower bound
        fill='toself',
        fillcolor='rgba(0,100,80,0.2)',  # Semi-transparent color
        line=dict(color='rgba(255,255,255,0)'),
        showlegend=False,
        name=f'{fig_confidence_interval}% Confidence Interval'
    ))

    # Calculate Spearman correlations
    rho1, p1 = stats.spearmanr([i for i,j in enumerate(x_values)], y_values)
    rho1, p1 = convert_r_p(rho1, p1)

    rho2, p2 = stats.spearmanr([i for i,j in enumerate(x_values)], loess_y)
    rho2, p2 = convert_r_p(rho2, p2)

    annotation_text = (
        f"<b>Spearman correlation</b><br>"
        f"Original data: ρ = {rho1} (p={p1})<br>"
        f"LOESS data: ρ = {rho2} (p={p2})"
    )

    fig.add_annotation(
        text=annotation_text,
        xref="paper",
        yref="paper",
        x=0,
        y=1,
        xanchor="left",
        yanchor="bottom",
        showarrow=False,
        align="left"
    )


    if fig_taxon == 'Hash':
        fig_taxon = st.session_state['TTT_hash'][:-1]

    # Default Layout
    fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize), title=f'{fig_taxon} Richness', rangemode='tozero')
    fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize), dtick='linear')
    fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)
    fig.update_xaxes(type="category",range=[-0.5, len(x_values) - 0.5])

    # Create df
    df = pd.DataFrame()
    df['Samples (in order)'] = x_values
    df[f'{fig_taxon} Richness'] = y_values
    df[f'{fig_taxon} Richness (LOESS)'] = loess_y

    # Plot
    st.plotly_chart(fig, config=st.session_state['TTT_config'])
    export_plot("Time_series", f"{name}_{fig_taxon}", fig, "plotly")
    export_table("Time_series", f"{name}_{fig_taxon}", df, "xlsx")

########################################################################################################################
# API
def gbif_accession():

    ## create copies of the dataframes
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    table_display_name = st.session_state['table_display_name']

    ## store results here
    species_dict = {'': ['', '', '', '', '']}
    OTU_species = []
    OTU_keys = []
    OTU_genusKeys = []
    OTU_synonyms = []
    OTU_iucnRedListCategory = []
    OTU_habitat = []
    OTU_gbif_link = []

    for row in stqdm(df_taxon_table.values.tolist()):
        species = row[7]

        ## request GBIF
        if species != '':
            if species not in species_dict.keys():
                query = species.replace(' ', '%20')
                # Initialize a counter for the number of attempts
                attempts = 0
                # Initialize a flag for whether the request was successful
                success = False

                while not success and attempts < 20:
                    try:
                        # Make a GET request to the API
                        response = requests.get(f'https://api.gbif.org/v1/species/match?name={query}')

                        # If the request was successful, parse the response text as JSON
                        if response.status_code == 200:
                            data = json.loads(response.text)
                            success = True
                        else:
                            attempts += 1
                    except requests.exceptions.RequestException as e:
                        # If there was a network problem (e.g. DNS resolution, refused connection, etc), increment the counter and try again
                        attempts += 1

                # If the request was not successful after 20 attempts, print a warning
                if not success:
                    species_dict[species] = [''] * 6
                else:
                    # add data to dict
                    species_dict[species] = [data.get('species', ''), data.get('speciesKey', ''),
                                             data.get('genusKey', ''), data.get('synonym', '')]

                    # Make a GET request to the API to get IUCN Red List Category
                    response = requests.get(
                        f"https://api.gbif.org/v1/species/{data.get('speciesKey', '')}/iucnRedListCategory")
                    if response.status_code == 200:
                        iucn_data = json.loads(response.text)
                        species_dict[species].append(iucn_data.get('code', ''))
                    else:
                        species_dict[species].append('')

                    # Make a GET request to the API to get habitat
                    response = requests.get(
                        f"https://api.gbif.org/v1/species/{data.get('speciesKey', '')}/speciesProfiles")
                    if response.status_code == 200:
                        habitat_data = json.loads(response.text)
                        habitat = ', '.join([j for j in sorted(set([i.get('habitat', '').lower() for i in habitat_data['results']])) if j != ''])
                        species_dict[species].append(habitat)
                    else:
                        species_dict[species].append('')

                time.sleep(0.5)

            # create link
            gbif_url = f'https://www.gbif.org/species/{str(species_dict[species][1])}'

            ## append results
            OTU_species.append(species_dict[species][0])
            OTU_keys.append(str(species_dict[species][1]))
            OTU_genusKeys.append(str(species_dict[species][2]))
            OTU_synonyms.append(str(species_dict[species][3]))
            OTU_iucnRedListCategory.append(species_dict[species][4])
            OTU_habitat.append(species_dict[species][5])
            OTU_gbif_link.append(gbif_url)
        else:
            ## append results
            OTU_species.append('')
            OTU_keys.append('')
            OTU_genusKeys.append('')
            OTU_synonyms.append('')
            OTU_iucnRedListCategory.append('')
            OTU_habitat.append('')
            OTU_gbif_link.append('')

    ## append to traits df
    similarity_loc = df_taxon_table.columns.tolist().index('Similarity')+1
    df_taxon_table.insert(loc=similarity_loc, column='speciesKey', value=OTU_keys)
    df_taxon_table.insert(loc=similarity_loc, column='genusKey', value=OTU_genusKeys)
    df_taxon_table.insert(loc=similarity_loc, column='GBIF species', value=OTU_species)
    df_taxon_table.insert(loc=similarity_loc, column='Synonym', value=OTU_synonyms)
    df_taxon_table.insert(loc=similarity_loc, column='iucnRedListCategory', value=OTU_iucnRedListCategory)
    df_taxon_table.insert(loc=similarity_loc, column='Habitat', value=OTU_habitat)
    df_taxon_table.insert(loc=similarity_loc, column='GBIF', value=OTU_gbif_link)

    ## export table
    export_taxon_table(table_display_name, df_taxon_table, df_metadata_table, '')
    st.session_state['df_taxon_table'] = df_taxon_table.copy()
    st.success('GBIF Data was implemented!')

def gbif_upload_conversion():
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    table_display_name = st.session_state['table_display_name']
    df_samples = st.session_state['df_samples']

    # User input
    gbif_target_gene = st.session_state['gbif_target_gene']
    gbif_fwd_primer_name = st.session_state['gbif_fwd_primer_name']
    gbif_fwd_primer_seq = st.session_state['gbif_fwd_primer_seq']
    gbif_rvs_primer_name = st.session_state['gbif_rvs_primer_name']
    gbif_rvs_primer_seq = st.session_state['gbif_rvs_primer_seq']

    # Metadata Input
    if gbif_date != 'None':
        gbif_date_col = st.session_state['gbif_date']
        gbif_date_values = df_metadata_table[gbif_date_col]
    if gbif_lat != 'None':
        gbif_lat_col = st.session_state['gbif_lat']
        gbif_lat_values = df_metadata_table[gbif_lat_col]
    if gbif_lon != 'None':
        gbif_lon_col = st.session_state['gbif_lon']
        gbif_lon_col_values = df_metadata_table[gbif_lon_col]

    # Create OTU Table (Sheet)
    OTU_table = df_taxon_table[df_samples]
    OTU_table.index = df_taxon_table['Hash']
    OTU_table.index.name = ''

    # Create Taxonomy Table (Sheet)
    taxonomy_table = df_taxon_table[['Hash', 'Seq', 'Species', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']].copy()
    taxonomy_table.columns = ['id', 'DNA_sequence', 'scientificName', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus']

    # Create Samples Table (Sheet)
    samples_table = pd.DataFrame(df_samples, columns=['id'])
    samples_table['decimalLatitude'] = gbif_lat_values
    samples_table['decimalLongitude'] = gbif_lon_col_values
    samples_table['eventDate'] = gbif_date_col

    # Create Study Table (Sheet)
    study_table = pd.DataFrame([['target_gene', gbif_target_gene],
                                ['pcr_primer_forward', gbif_fwd_primer_seq],
                                ['pcr_primer_name_forward', gbif_fwd_primer_name],
                                ['pcr_primer_reverse', gbif_rvs_primer_seq],
                                ['pcr_primer_name_reverse', gbif_rvs_primer_name],
                                ], columns=['term', 'value'])

    active_project_path = st.session_state['active_project_path']
    folder_name = 'GBIF'
    os.makedirs(active_project_path / folder_name, exist_ok=True)
    export_file = active_project_path / folder_name / f"GBIF_Upload_{table_display_name.name}"
    with pd.ExcelWriter(export_file, engine="xlsxwriter") as writer:
        OTU_table.to_excel(writer, sheet_name='OTU_table', index=True)
        taxonomy_table.to_excel(writer, sheet_name='Taxonomy', index=False)
        samples_table.to_excel(writer, sheet_name='Samples', index=False)
        study_table.to_excel(writer, sheet_name='Study', index=False)

    st.success('GBIF Upload Table was created!')
    open_file(export_file)

def gbif_taxonomy_validation():
    pass

def gbif_occurrence_download():
    df_taxon_table = st.session_state['df_taxon_table'].copy()

    username = 'till_macher'
    password = 'ZN_kw7@!gBsABNn'
    email = 'macher@uni-trier.de'

    occ.download_get(key='0024088-260226173443078')

    active_project_path = st.session_state['active_project_path']
    folder_name = 'GBIF'
    os.makedirs(active_project_path / folder_name, exist_ok=True)
    specieskey_json = active_project_path / folder_name / f"GBIF_speciesKeys.json"
    records_dict = {}
    if specieskey_json.exists():
        with open(specieskey_json, "r") as f:
            records_dict = json.load(f)

    if 'speciesKey' in df_taxon_table.columns.tolist():
        unique_keys = [i for i in df_taxon_table['speciesKey'].unique() if i != '']
        for specieskey in unique_keys:
            specieskey = int(specieskey)
            if str(specieskey) not in records_dict.keys():
                try:
                    res = occ.search(taxonKey=specieskey, limit=300)
                    records_dict[str(specieskey)] = res['results']
                except:

                    records_dict[str(specieskey)] = ''
    else:
        st.warning('Please run the GBIF basics tool first!')

    with open(specieskey_json, "w") as f:
        json.dump(records_dict, f)

########################################################################################################################
# Map Modules
def sample_map():

    # Style
    fontsize = st.session_state['TTT_fontsize']
    template = st.session_state['TTT_template']
    showlegend = st.session_state['TTT_showlegend']
    scattersize = st.session_state['TTT_scattersize']
    showlabels = st.session_state['nmds_showlabels']
    linewidth = st.session_state['TTT_linewidth']
    TTT_color1 = st.session_state['TTT_color1']
    TTT_color2 = st.session_state['TTT_color2']
    map_style = st.session_state['map_style']
    sample_radius = st.session_state['sample_radius']
    map_zoom = st.session_state['map_zoom']
    map_groups = st.session_state['map_groups']
    fig_color = st.session_state['map_color']
    map_save = st.session_state['map_save']

    name = st.session_state['table_display_name'].stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples']
    table_display_name = st.session_state['table_display_name']

    # Colors
    groups = df_metadata_table[map_groups].unique()
    n_groups = len(groups)
    if fig_color == 'Color Scale':
        colors = get_colors_from_scale(TTT_colorscale, n_groups)
    if fig_color == 'Color Sequence':
        colors = get_colors_from_sequence(TTT_colorsequence, n_groups)
    if fig_color == 'Single Color':
        colors = [TTT_color1] * n_groups
    df_metadata_table = colors_to_metadata(df_metadata_table, map_groups, colors)

    ####################################################################################################
    # Prepare site table
    lon_col = "Longitude"
    lat_col = "Latitude"
    df_sites = df_metadata_table[["Samples", lat_col, lon_col]].copy()
    df_sites.columns = ["site", "lat", "lon"]
    # Drop samples without Lat Long
    df_sites = df_sites[(df_sites['lat'] != '') & (df_sites['lon'] != '')]
    df_sites['color'] = df_metadata_table['TTT_metadatacolor']
    if len(df_sites) == 0:
        return

    ####################################################################################################
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(df_sites, geometry=gpd.points_from_xy(df_sites.lon, df_sites.lat), crs="EPSG:4326")
    ####################################################################################################
    # Project to meters for buffering
    gdf_m = gdf.to_crs("EPSG:3035")
    # Create circular buffers
    gdf_m["buffer"] = gdf_m.buffer(sample_radius * 1000)
    ####################################################################################################
    # Convert buffers back to WGS84
    buffers_wgs = gdf_m.set_geometry("buffer").to_crs("EPSG:4326")
    ####################################################################################################
    # Merge overlapping buffers
    merged_buffers = buffers_wgs["buffer"].union_all()
    # Optional simplification for faster plotting
    merged_buffers = merged_buffers.simplify(0.001)
    ####################################################################################################
    # Handle Polygon / MultiPolygon
    if isinstance(merged_buffers, Polygon):
        polygons = [merged_buffers]
    elif isinstance(merged_buffers, MultiPolygon):
        polygons = list(merged_buffers.geoms)
    else:
        polygons = []
    ####################################################################################################
    # Create plot
    fig = go.Figure()
    # Add merged buffer polygons
    showlegend = False
    for poly in polygons:
        lon_outline, lat_outline = poly.exterior.coords.xy
        fig.add_trace(go.Scattermap(
                lat=list(lat_outline),
                lon=list(lon_outline),
                mode="lines", fill="toself", fillcolor="rgba(0,0,0,0.1)", line=dict(width=2, color=TTT_color2),
                name="Sampling buffer", showlegend=showlegend))
    ####################################################################################################
    # Add sample points
    fig.add_trace(
        go.Scattermap(
            lon=df_sites["lon"],
            lat=df_sites["lat"],
            text=df_sites["site"],
            mode="markers+text",
            marker=dict(size=scattersize*1.1, color=df_sites["color"]),
            textfont=dict(size=fontsize, color="black"),
            textposition="top center", showlegend=False,
            name="Samples"))
    ####################################################################################################
    # Map layout
    fig.update_layout(
        map=dict(
            style=map_style,
            zoom=map_zoom,
            center=dict(lat=df_sites["lat"].mean(), lon=df_sites["lon"].mean())),
            title=f"{table_display_name.stem} – Sampling Sites ({sample_radius} km)")
    ####################################################################################################
    # Show map in Streamlit
    st.plotly_chart(fig, width='stretch', height=TTT_height, config=st.session_state['TTT_config'])

    if map_save == 'Yes':
        export_plot("Maps", f"{name}_map", fig, "plotly")

def gbif_occurrence_validation():

    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples']
    table_display_name = st.session_state['table_display_name']
    name = table_display_name.stem
    sample_radius = st.session_state['sample_radius']  # km

    # check if table exists
    active_project_path = st.session_state['active_project_path']
    folder_name = 'Maps'
    suffix = f"{name}_gbif_{sample_radius}km"
    file_xlsx = active_project_path / folder_name / f'{suffix}.xlsx'

    if file_xlsx.exists():
        st.info('Occurrence Analysis was already conducted: Please delete the output file to rerun the analysis!')
        return

    # Prepare sites
    lon_col, lat_col = "Longitude", "Latitude"
    df_sites = df_metadata_table[["Samples", lat_col, lon_col]].copy()
    df_sites.columns = ["site", "lat", "lon"]
    df_sites = df_sites[(df_sites['lat'] != '') & (df_sites['lon'] != '')]

    # Convert to GeoDataFrame
    gdf = gpd.GeoDataFrame(df_sites, geometry=gpd.points_from_xy(df_sites.lon, df_sites.lat), crs="EPSG:4326")

    # Project to meters, create circular buffers
    gdf_m = gdf.to_crs("EPSG:3035")
    gdf_m["buffer"] = gdf_m.buffer(sample_radius * 1000)

    # Convert back to WGS84, simplify, and force counter-clockwise
    def buffer_to_wkt(buf):
        buf_simplified = buf.simplify(0.01, preserve_topology=True)
        buf_ccw = orient(buf_simplified, sign=1.0)
        return buf_ccw.wkt

    buffers_wgs = gdf_m.set_geometry("buffer").to_crs("EPSG:4326")
    df_sites["polygon"] = buffers_wgs["buffer"].apply(buffer_to_wkt)

    # Initialize occurrence counts per species
    species_keys = [i for i in df_taxon_table['Species'].unique() if i != '']
    res = []
    for species_name in stqdm(species_keys, desc='GBIF occurrence query'):
        counts = [species_name]
        for _, row in df_sites.iterrows():
            try:
                gbif_res = occ.search(scientificName=species_name, geometry=row["polygon"], limit=0)
                counts.append(gbif_res["count"])
            except Exception:
                counts.append(0)
            time.sleep(0.1)
        res.append(counts)
    res_df = pd.DataFrame(res, columns=['Species'] + df_samples)

    # calculate correlation to identify species which are unlikely to occur
    df_taxon_table_s = simple_taxon_table(df_taxon_table)[['Species'] + df_samples]
    df_taxon_table_s = df_taxon_table_s[df_taxon_table_s['Species'] != '']
    df_relative = df_taxon_table_s.copy()
    df_relative[df_samples] = df_relative[df_samples].div(df_relative[df_samples].sum(axis=0), axis=1)
    df_relative[df_samples] = df_relative[df_samples] * 100
    df_relative = df_relative.round(4)

    # Melt both tables to long format
    gbif_long = res_df.melt(id_vars='Species', var_name='Sample', value_name='GBIF (occ.)')
    reads_long = df_relative.melt(id_vars='Species', var_name='Sample', value_name='Reads (%)')
    # Merge by Species and Sample
    combined_long = pd.merge(reads_long, gbif_long, on=['Species', 'Sample'])
    # Optional: reorder columns
    combined_long = combined_long[['Sample', 'Species', 'Reads (%)', 'GBIF (occ.)']]
    combined_long = combined_long.sort_values(['Species', 'Sample'])

    # Export Table
    export_table("Maps", suffix, combined_long, "xlsx", False)

def taxon_distribution_map():
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples']
    table_display_name = st.session_state['table_display_name']
    name = table_display_name.stem
    distmap_level = st.session_state['distmap_level']
    distmap_taxon = st.session_state['distmap_taxon']
    TTT_colorscale = st.session_state['TTT_colorscale']
    map_style = st.session_state['distmap_style']
    map_zoom = st.session_state['distmap_zoom']
    distmap_save = st.session_state['distmap_save']
    distmap_measure = st.session_state['distmap_measure']
    TTT_scattersize = st.session_state['TTT_scattersize']

    ####################################################################################################
    # Filter to taxon to display
    if distmap_measure == 'Rel. Reads':
        df_taxon_table_s = simple_taxon_table(df_taxon_table)
        df_relative = df_taxon_table_s.copy()
        df_relative[df_samples] = df_relative[df_samples].div(df_relative[df_samples].sum(axis=0), axis=1)
        df_relative[df_samples] = df_relative[df_samples] * 100
        alpha_div_values = df_relative[df_relative[distmap_level] == distmap_taxon][df_samples].sum().values.tolist()
        title=f'{distmap_taxon} (Rel. Reads)'
        cmin, cmax = 0, 100
    else:
        df_taxon_table_s = df_taxon_table[df_taxon_table[distmap_level] == distmap_taxon]
        alpha_div_values = []
        for sample in df_samples:
            n_hashes = len(df_taxon_table_s[df_taxon_table_s[sample] != 0])
            alpha_div_values.append(n_hashes)
        title=f'{distmap_taxon} ({st.session_state["TTT_hash"]})'
        cmin, cmax = 1, max(alpha_div_values)+1

    ####################################################################################################
    # Prepare site table
    lon_col = "Longitude"
    lat_col = "Latitude"
    df_sites = df_metadata_table[["Samples", lat_col, lon_col]].copy()
    df_sites.columns = ["site", "lat", "lon"]
    # Drop samples without Lat Long
    df_sites = df_sites[(df_sites['lat'] != '') & (df_sites['lon'] != '')]
    df_sites[distmap_taxon] = alpha_div_values
    if len(df_sites) == 0:
        return

    ####################################################################################################
    # Create figure
    fig = go.Figure()

    # add no detection samples
    fig_df = df_sites[df_sites[distmap_taxon] == 0.0]
    fig.add_trace(go.Scattermap(
        lat=fig_df['lat'],
        lon=fig_df['lon'],
        mode='markers+text',
        text=fig_df['site'],
        textposition='top center',
        showlegend=False,
        hovertext=fig_df[distmap_taxon],
        marker=dict(
            size=TTT_scattersize*2,  # marker size
            color='white',  # the column you want to color by
            opacity=0.3,
        )
    ))

    # add detection samples
    fig_df = df_sites[df_sites[distmap_taxon] > 0]
    fig.add_trace(go.Scattermap(
        lat=fig_df['lat'],
        lon=fig_df['lon'],
        mode='markers+text',
        text=fig_df['site'],
        textposition='top center',
        line=dict(width=1, color='black'),
        showlegend=False,
        hovertext=fig_df[distmap_taxon],
        marker=dict(
            size=TTT_scattersize*2,  # marker size
            color=fig_df[distmap_taxon],
            colorscale=TTT_colorscale,
            colorbar=dict(title=title),
            cmin=cmin,
            cmax=cmax,
            opacity=0.9
        )
    ))

    ####################################################################################################
    # Map layout
    fig.update_layout(
        map=dict(
            style=map_style,
            zoom=map_zoom,
            center=dict(lat=df_sites["lat"].mean(), lon=df_sites["lon"].mean())),
            title=f"{table_display_name.stem} – {distmap_taxon}")
    ####################################################################################################
    st.plotly_chart(fig, width='stretch', height=TTT_height, config=st.session_state['TTT_config'])

    if distmap_save == 'Yes':
        export_plot("Maps", f"{name}_{distmap_taxon}_{distmap_measure}_distmap", fig, "plotly")

########################################################################################################################
# ESC calculation

def convert_to_perlodes(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, tool_settings):
    # Make copies of input dataframes to prevent altering the originals
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    river_types_dict = {i:j for i,j in metadata_df[['Sample', 'Perlodes_river_type']].values.tolist()}
    taxa_list_dict = {i:j for i,j in metadata_df[['Sample', 'Perlodes_taxa_list']].values.tolist()}
    usage_dict = {i:j for i,j in metadata_df[['Sample', 'Perlodes_usage']].values.tolist()}

    # Ensure all samples have metadata in all three categories
    samples = [
        i for i in samples
        if i in river_types_dict and pd.notna(river_types_dict[i]) and
           i in taxa_list_dict and pd.notna(taxa_list_dict[i]) and
           i in usage_dict and pd.notna(usage_dict[i])
            ]

    if len(samples) == 0:
        st.warning('Please fill out the required metadata for at least one sample!')
        return

    # Extract relevant settings from the tool settings
    presence_absence = tool_settings['presence_absence']
    taxon_table_taxonomy = taxon_table_df.columns.tolist()[1:7]

    all_reads = []
    for sample in samples:
        total_reads = taxon_table_df[sample].sum()
        taxon_dict = {}
        for OTU in taxon_table_df[taxon_table_taxonomy + [sample]].values.tolist():
            taxon = [item for item in OTU[0:6] if item != ''][-1]
            if taxon not in taxon_dict.keys():
                n_reads = OTU[-1]
                taxon_dict[taxon] = n_reads
            else:
                n_reads = OTU[-1]
                taxon_dict[taxon] = taxon_dict[taxon] + n_reads
        if presence_absence == True:
            pa_list = [1 if i != 0 else 0 for i in taxon_dict.values()]
            all_reads.append(pa_list)
        elif presence_absence == 'Relative':
            rel_list = [round(i / total_reads * 100, 2) for i in taxon_dict.values()]
            all_reads.append(rel_list)
        else:
            all_reads.append(list(taxon_dict.values()))

    # create initial df
    index = [i for i in range(1,1+len(taxon_dict.keys()))]
    species = list(taxon_dict.keys())
    perlodes_df = pd.DataFrame(species, columns=['species'])
    perlodes_df.insert(0, '', index)
    perlodes_df = pd.concat([perlodes_df, pd.DataFrame(all_reads, index=samples).transpose()], axis=1)

    # convert the df to match the (complicated) perlodes input format
    columns = ['ID_ART', 'TAXON_NAME'] + samples
    gewässertyp = ['Gewässertyp', ''] + [river_types_dict[i] for i in samples]
    Taxaliste = ['Taxaliste', ''] + [taxa_list_dict[i] for i in samples]
    Nutzung = ['Nutzung', ''] + [usage_dict[i] for i in samples]
    taxa = [[1] + i[1:] for i in perlodes_df.values.tolist()]
    df_list = [gewässertyp] + [Taxaliste] + [Nutzung] + taxa
    perlodes_input_df = pd.DataFrame(df_list, columns=columns)

    # write the filtered list to a dataframe
    if presence_absence == True:
        perlodes_directory = Path(str(path_to_outdirs) + "/" + "Perlodes" + "/" + taxon_table_xlsx.stem)
        perlodes_xlsx = Path(str(perlodes_directory) + "_Perlodes_PA.xlsx")
        perlodes_input_df.to_excel(perlodes_xlsx, index=False)
    elif presence_absence == 'Relative':
        perlodes_directory = Path(str(path_to_outdirs) + "/" + "Perlodes" + "/" + taxon_table_xlsx.stem)
        perlodes_xlsx = Path(str(perlodes_directory) + "_Perlodes_RELATIVE_ABUNDANCE.xlsx")
        perlodes_input_df.to_excel(perlodes_xlsx, index=False)
    else:
        perlodes_directory = Path(str(path_to_outdirs) + "/" + "Perlodes" + "/" + taxon_table_xlsx.stem)
        perlodes_xlsx = Path(str(perlodes_directory) + "_Perlodes_ABUNDANCE.xlsx")
        perlodes_input_df.to_excel(perlodes_xlsx, index=False)

    st.success(f'Wrote Perlodes input file to: {perlodes_xlsx}')

def convert_to_phylib(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, tool_settings):

    #get the taxonomy from the operational taxon list
    operational_taxon_list_df = pd.read_excel(Path(operational_taxon_list), sheet_name="TTT import").fillna('')

    # load the taxon table and create a list
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    meta_data_df = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)
    TaXon_table_taxonomy = TaXon_table_df.columns.tolist()[0:7]
    samples_list = TaXon_table_df.columns.tolist()[10:]

    ## load the metadata -> freshwater type
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
    metadata_loc = Meta_data_table_df.columns.tolist().index(meta_data_to_test)
    types_dict = {i[0]:i[1] for i in Meta_data_table_df[['Samples', meta_data_to_test]].values.tolist()}

    ## drop samples with metadata called nan (= empty)
    drop_samples = [i[0] for i in Meta_data_table_df.values.tolist() if i[metadata_loc] == ""]

    ## test if samples have metadata
    if drop_samples != []:
        sg.PopupError("Please fill out all the metadata for all samples.")

    ## test if all metadata for Phylib is available
    elif len(set([True if i in Meta_data_table_df.columns.tolist() else False for i in ['Ökoregion', 'Makrophytenverödung', 'Begründung', 'Helophytendominanz', 'Diatomeentyp', 'Phytobenthostyp', 'Makrophytentyp', 'WRRL-Typ', 'Gesamtdeckungsgrad']])) != 1:
        sg.PopupError("Please fill out all the required phylib metadata for all samples.")

    elif sorted(Meta_data_table_samples) == sorted(samples_list):
        # store hits and dropped OTUs
        hit_list, dropped_list = [], []

        # loop through the taxon table
        for taxonomy in TaXon_table_df[TaXon_table_taxonomy].values.tolist():
            ## test species
            if taxonomy[6] != '' and taxonomy[6] in operational_taxon_list_df['Species'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Species'].str.contains(taxonomy[6])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test genus
            elif taxonomy[5] != '' and taxonomy[5] in operational_taxon_list_df['Genus'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Genus'].str.contains(taxonomy[5])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test family
            elif taxonomy[4] != '' and taxonomy[4] in operational_taxon_list_df['Family'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Family'].str.contains(taxonomy[4])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test order
            elif taxonomy[3] != '' and taxonomy[3] in operational_taxon_list_df['Order'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Order'].str.contains(taxonomy[3])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test class
            elif taxonomy[2] != '' and taxonomy[2] in operational_taxon_list_df['Class'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Class'].str.contains(taxonomy[2])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test phylum
            elif taxonomy[1] != '' and taxonomy[1] in operational_taxon_list_df['Phylum'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Phylum'].str.contains(taxonomy[1])][['DV-NR.', 'Taxon']].values.tolist()[0]
            else:
                res = taxonomy + ['', '']
                dropped_list.append(taxonomy)

            hit_list.append(res)

        ## create a dataframe
        ## create a new df and export it as TaXon table
        df = pd.DataFrame(hit_list, columns=TaXon_table_df.columns[0:7].values.tolist() + ['DV-NR.', 'Taxon (Phylib)'])
        concatenated_df = pd.concat([df, TaXon_table_df[samples_list + ['Similarity', 'Status', 'seq']]], axis=1)
        reordered_df = concatenated_df[['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'Status', 'DV-NR.', 'Taxon (Phylib)', 'seq'] + samples_list]
        reordered_metadata_df = add_metadata(reordered_df, meta_data_df)
        phylib_TaXon_Table = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem + "_phylib.xlsx")
        reordered_metadata_df.to_excel(phylib_TaXon_Table, index=False)

        # create an output list for phylib
        messwerte_list = []
        messtellen_list = []

        for sample in samples_list:
            ## collect information about the sampling site
            sample_metadata_df = Meta_data_table_df.loc[Meta_data_table_df['Samples'] == sample]
            messstelle = sample
            ökoregion = sample_metadata_df['Ökoregion'].values.tolist()[0]
            Makrophytenverödung = sample_metadata_df['Makrophytenverödung'].values.tolist()[0]
            Begründung = sample_metadata_df['Begründung'].values.tolist()[0]
            Helophytendominanz = sample_metadata_df['Helophytendominanz'].values.tolist()[0]
            Diatomeentyp = sample_metadata_df['Diatomeentyp'].values.tolist()[0]
            Phytobenthostyp = sample_metadata_df['Phytobenthostyp'].values.tolist()[0]
            Makrophytentyp = sample_metadata_df['Makrophytentyp'].values.tolist()[0]
            WRRL_Typ = sample_metadata_df['WRRL-Typ'].values.tolist()[0]
            Gesamtdeckungsgrad = sample_metadata_df['Gesamtdeckungsgrad'].values.tolist()[0]
            messtellen_list.append([messstelle, ökoregion, Makrophytenverödung, Begründung, Helophytendominanz, Diatomeentyp, Phytobenthostyp, Makrophytentyp, WRRL_Typ, Gesamtdeckungsgrad])

            ## calcualte sum of reads for taxa with multiple OTUs
            sample_df = reordered_metadata_df[[sample, 'DV-NR.', 'Taxon (Phylib)']]
            reads_dict = {}
            for taxon in sample_df.values.tolist():
                key = taxon[2]
                values = taxon[:2]
                if key != '':
                    if key not in reads_dict.keys():
                        reads_dict[key] = values
                    else:
                        reads_dict[key] = [reads_dict[key][0] + values[0], values[1]]

            ## remove duplicates
            samples_taxa = [[key]+values for key,values in reads_dict.items() if values[0] != 0]

            ## if PA data is required: Convert to 1/0
            if presence_absence == True:
                samples_taxa = [[i[0], 1, i[2]] for i in samples_taxa]

            ## calculate the overall number of reads/specimens. This is required for later relative abundance calculation
            sum_measurement = sum([i[1] for i in samples_taxa])

            ## loop through all taxa
            for taxon in samples_taxa:
                ## if the taxon (assigned to the OTU) is present in the sample and present on the OTL, continue
                ## add all relevant information to list (for df)
                probe = sample
                taxon_id = taxon[2]
                taxonname = taxon[0]
                form = "o.A."
                messwert = taxon[1] / sum_measurement * 100
                einheit = "%"
                cf = ""
                messwerte_list.append([messstelle, probe, taxon_id, taxonname, form, messwert, einheit, cf])

        phylib_df_1 = pd.DataFrame(messtellen_list, columns=["Messstelle", "Ökoregion", "Makrophytenverödung", "Begründung", "Helophytendominanz", "Diatomeentyp", "Phytobenthostyp", "Makrophytentyp", "WRRL-Typ", "Gesamtdeckungsgrad"])
        phylib_df_2 = pd.DataFrame(messwerte_list, columns=["Messstelle", "Probe", "Taxon", "Taxonname", "Form", "Messwert", "Einheit", "cf"])

        if presence_absence == False:
            phylib_directory = Path(str(path_to_outdirs) + "/" + "Phylib" + "/" + TaXon_table_xlsx.stem)
            phylib_xlsx = Path(str(phylib_directory) + "_phylib_ABUNDANCE.xlsx")
            writer = pd.ExcelWriter(phylib_xlsx, engine='xlsxwriter')
        else:
            phylib_directory = Path(str(path_to_outdirs) + "/" + "Phylib" + "/" + TaXon_table_xlsx.stem)
            phylib_xlsx = Path(str(phylib_directory) + "_phylib_PA.xlsx")
            writer = pd.ExcelWriter(phylib_xlsx, engine='xlsxwriter')
        phylib_df_1.to_excel(writer, sheet_name='Messstelle', index=False)
        phylib_df_2.to_excel(writer, sheet_name='Messwerte', index=False)
        writer.save()

        ################################################################################################################
        ## create some plots and provide statistics
        unique_original_taxa = list(k for k, _ in itertools.groupby(reordered_metadata_df[['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].values.tolist()))
        len(unique_original_taxa)

        unique_OTL_taxa = [i for i in set(reordered_metadata_df['Taxon (Phylib)'].values.tolist()) if i != '']
        len(unique_OTL_taxa)

        phylib_conversion_loss = {}
        for taxon in unique_OTL_taxa:
            all_taxa = reordered_metadata_df[reordered_metadata_df['Taxon (Phylib)'].str.contains(taxon)][['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].values.tolist()
            unique_taxa = set([''.join(i) for i in all_taxa])
            phylib_conversion_loss[taxon] = len(unique_taxa)

        perlodes_conversion_loss = {k: v for k, v in sorted(phylib_conversion_loss.items(), key=lambda item: item[1])}

        fig = go.Figure()
        x_values = [i[:19] for i in list(perlodes_conversion_loss.keys())[-30:][::-1]] ## [:20] to cut names off
        y_values = list(perlodes_conversion_loss.values())[-30:][::-1]
        fig.add_trace(go.Bar(x=x_values, y=y_values, marker_color='blue'))
        fig.update_yaxes(title='# taxa')
        fig.update_layout(title='', showlegend=False, font_size=14, width=width_value, height=height_value, template=template)
        phylib_plot = Path(str(phylib_directory) + "_phylib.pdf")
        fig.write_image(phylib_plot)
        phylib_plot = Path(str(phylib_directory) + "_phylib.html")
        fig.write_html(phylib_plot)

def calculate_diathor_index():
    # Tables & Samples
    table_display_name = st.session_state['table_display_name']
    name = table_display_name.stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples']

    # Variables
    diathor_measure = st.session_state['diathor_measure']

    # collect species
    # Remove empty species once
    df = df_taxon_table[df_taxon_table['Species'] != ''].copy()

    # Precompute total reads per sample (only once!)
    total_reads = df[df_samples].sum()

    # Aggregate reads per species
    grouped = df.groupby('Species')[df_samples].sum()

    if diathor_measure == 'Presence/Absence':
        result = (grouped > 0).astype(int)

    elif diathor_measure == 'Rel. Reads':
        result = grouped.div(total_reads, axis=1) * 100

    else:  # hash count
        result = (
            df.groupby('Species')[df_samples]
            .apply(lambda x: (x != 0).sum())
        )

    # Sort species
    result = result.sort_index()

    # Add index column like before
    result.insert(0, 'species', result.index)

    # Reset index for clean df
    diathor_df = result.reset_index(drop=True)

    # drop empty samples
    sample_sums = diathor_df[df_samples].sum()
    samples_to_keep = sample_sums[sample_sums != 0].index.tolist()
    diathor_df = diathor_df[['species'] + samples_to_keep]

    # write the filtered list to a dataframe
    table_name = f"{diathor_measure.replace('/', '_').replace(' ', '')}"
    diathor_csv = export_table("Diathor", table_name, diathor_df, "csv")
    st.success(f'Wrote Diathor input file.')

    ## Calculate Diathor
    st.info('Diathor output progress is printed to the terminal.')
    output_folder = diathor_csv.parent.joinpath(diathor_csv.stem)
    os.makedirs(output_folder, exist_ok=True)
    output_xlsx = output_folder.joinpath(str(diathor_csv.name).replace('.csv', '.xlsx'))
    log_txt = Path(str(output_xlsx).replace('.xlsx', '.log'))
    f = open(log_txt, 'w')
    path_to_ttt = Path(__file__).resolve().parent
    script_path = path_to_ttt / 'Rscripts' / 'diathor.R'
    subprocess.call(['Rscript', script_path, '-i', diathor_csv, '-o', output_xlsx], stderr=f)
    f.close()
    st.success(f'Finished Diathor calculations.')

def calculate_efi_index():

    # Tables & Samples
    table_display_name = st.session_state['table_display_name']
    name = table_display_name.stem
    df_taxon_table = st.session_state['df_taxon_table'].copy()
    df_metadata_table = st.session_state['df_metadata_table'].copy()

    # Variables
    efi_placeholder = st.session_state['efi_placeholder']
    efi_measure = st.session_state['efi_measure']
    taxon_table_taxonomy = list(st.session_state['tax_level_plurals'].keys())

    # Check if all required columns are present in the table:
    efi_columns = {
        "EFI_Day": 1,
        "EFI_Month": 1,
        "EFI_Year": 2026,
        "EFI_Longitude": 39.664704931322646,
        "EFI_Latitude": -5.939361568912182,
        "EFI_Actual.river.slope": 1.1,
        "EFI_Temp.jul": 20,
        "EFI_Temp.jan": 10,
        "EFI_Floodplain.site": "No",
        "EFI_Water.source.type": "Pluvial",
        "EFI_Geomorph.river.type": "Meand regular",
        "EFI_Distance.from.source": 58,
        "EFI_Area.ctch": 450,
        "EFI_Natural.sediment": "Gravel/Pebble/Cobble",
        "EFI_Ecoreg": "Central highlands",
        "EFI_Eft.type": "T.15",
        "EFI_Fished.area": 500,
        "EFI_Method": "Wading"
    }
    available_metadata = df_metadata_table.columns.tolist()
    if not all(col in available_metadata for col in efi_columns.keys()):
        st.warning('Could not find all required EFI columns in the Metadata Table.')
        for col, values in efi_columns.items():
            if efi_placeholder == 'Yes':
                df_metadata_table[col] = values
            else:
                df_metadata_table[col] = ''
        export_taxon_table(table_display_name, df_taxon_table, df_metadata_table, '')
        st.info('You table was updated - please fill out and reload you Taxon Table.')
        open_file(table_display_name)
        return

    df_samples = st.session_state['df_samples'].copy()
    df_metadata_table_efi = df_metadata_table[['Samples'] + list(efi_columns.keys())]

    # Filter samples
    samples = [sample for sample in df_samples if '' not in df_metadata_table_efi[df_metadata_table_efi['Samples'] == sample].values.tolist()[0]]

    if len(samples) == 0:
        st.warning('Please fill out the required metadata for at least one sample!')
        return

    # efi_measure = 'Diversity'
    # collect species
    i = 1
    efi_df_values = []
    for sample in df_samples:
        if efi_measure == 'Presence/Absence' or efi_measure == 'Rel. Reads':
            taxon_dict = {}
            total_reads = df_taxon_table[sample].sum()
            for OTU in df_taxon_table[taxon_table_taxonomy + [sample]].values.tolist():
                taxon = [item for item in OTU[0:7] if item != ''][-1]
                if taxon not in taxon_dict.keys():
                    n_reads = OTU[-1]
                    taxon_dict[taxon] = n_reads
                else:
                    n_reads = OTU[-1]
                    taxon_dict[taxon] = taxon_dict[taxon] + n_reads
        else:
            taxa = [i for i in df_taxon_table[df_taxon_table[sample] != 0]['Species'].values.tolist() if i != '']
            taxon_dict = {i: taxa.count(i) for i in set(taxa)}

        efi_sample_metadata = df_metadata_table_efi[df_metadata_table_efi['Samples'] == sample].values.tolist()
        for taxon, reads in taxon_dict.items():
            if reads != 0:
                if efi_measure == 'Presence/Absence':
                    total_number_run1 = 1
                    number_length_below_150 = 1
                    number_length_above_150 = 0
                elif efi_measure == 'Rel. Reads':
                    total_number_run1 = round(reads / total_reads * 100, 3)
                    number_length_below_150 = total_number_run1
                    number_length_above_150 = 0
                else:
                    total_number_run1 = reads
                    number_length_below_150 = reads
                    number_length_above_150 = 0
                species_values = efi_sample_metadata[0] + [0, taxon, total_number_run1, number_length_below_150, number_length_above_150, 'abc', i, taxon]
                efi_df_values.append(species_values)
                i +=1

    # construct EFI dataframe
    efi_df = pd.DataFrame(efi_df_values, columns=['Sample.code'] + [i.replace('EFI_', '') for i in list(efi_columns.keys())] + ['Medit', 'Species', 'Total.number.run1', 'Number.length.below.150', 'Number.length.over.150', 'Sampling.location', 'code', 'species'])

    # add dummy values
    # other EFI will crash in some cases...
    path_to_ttt = Path(__file__).resolve().parent
    efi_table = path_to_ttt / 'WFD_conversion' / 'efi_table.xlsx'
    efi_dummy_df = pd.read_excel(efi_table, sheet_name='Sheet3').fillna('')
    efi_dummy_df = efi_dummy_df[efi_df.columns.tolist()]
    efi_df = pd.concat([efi_df, efi_dummy_df], ignore_index=True)
    efi_df['code'] = [i for i in range(1, len(efi_df)+1)]

    # write the filtered list to a dataframe
    table_name = f"{efi_measure.replace('/', '_').replace(' ', '')}"
    efi_xlsx = export_table("EFI", table_name, efi_df, "xlsx")
    st.success(f'Wrote EFI input file.')

    ####################################################################################################################
    ## Calculate EFI
    st.info('EFI output progress is printed to the terminal.')
    active_project_path = st.session_state['active_project_path']
    os.makedirs(active_project_path / "EFI", exist_ok=True)
    efi_results_xlsx = active_project_path / "EFI" / f"{name}_{table_name}_results.xlsx"
    efi_log = active_project_path / "EFI" / f"{name}_{table_name}.log"
    f = open(efi_log, 'w')
    script_path = path_to_ttt / 'Rscripts' / 'EFI.R'
    efi_directory = path_to_ttt / 'Rscripts'
    subprocess.call([
        'Rscript',
        script_path,
        '-d', efi_directory,
        '-i', efi_xlsx,
        '-o', efi_results_xlsx
    ], stderr=f)
    f.close()
    st.success(f'Finished EFI calculations!')

########################################################################################################################
# side bar

def TTT_variables():

    directories_to_create = ["Venn_diagrams","TaXon_tables", "Rarefaction_curves",
                             "Site_occupancy_plots", "Read_proportions_plots",
                             "Krona_charts", "Alpha_diversity", "Beta_diversity",
                             "PCoA_plots", "Replicate_analysis", "GBIF", "Occurrence_analysis",
                             "Per_taxon_statistics", "NMDS_plots", "Table_comparison", "Perlodes",
                             "EFI", "Diathor", "Phylib", "Time_series", "Basic_stats", "Fasta", "Import",
                             "Rarefaction_curves"]

    available_templates_list = ['TaxonTableTools', 'seaborn', 'ggplot2', 'simple_white', 'plotly', 'plotly_dark', 'presentation', 'plotly_white']

    available_clustering_units = ['OTUs', 'zOTUs', 'ESVs', 'ASVs']

    plotly_colors = ["aliceblue", "antiquewhite", "aqua", "aquamarine", "azure",
    "beige", "bisque", "black", "blanchedalmond", "blue",
    "blueviolet", "brown", "burlywood", "cadetblue",
    "chartreuse", "chocolate", "coral", "cornflowerblue",
    "cornsilk", "crimson", "cyan", "darkblue", "darkcyan",
    "darkgoldenrod", "darkgray", "darkgrey", "darkgreen",
    "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange",
    "darkorchid", "darkred", "darksalmon", "darkseagreen",
    "darkslateblue", "darkslategray", "darkslategrey",
    "darkturquoise", "darkviolet", "deeppink", "deepskyblue",
    "dimgray", "dimgrey", "dodgerblue", "firebrick",
    "floralwhite", "forestgreen", "fuchsia", "gainsboro",
    "ghostwhite", "gold", "goldenrod", "gray", "grey", "green",
    "greenyellow", "honeydew", "hotpink", "indianred", "indigo",
    "ivory", "khaki", "lavender", "lavenderblush", "lawngreen",
    "lemonchiffon", "lightblue", "lightcoral", "lightcyan",
    "lightgoldenrodyellow", "lightgray", "lightgrey",
    "lightgreen", "lightpink", "lightsalmon", "lightseagreen",
    "lightskyblue", "lightslategray", "lightslategrey",
    "lightsteelblue", "lightyellow", "lime", "limegreen",
    "linen", "magenta", "maroon", "mediumaquamarine",
    "mediumblue", "mediumorchid", "mediumpurple",
    "mediumseagreen", "mediumslateblue", "mediumspringgreen",
    "mediumturquoise", "mediumvioletred", "midnightblue",
    "mintcream", "mistyrose", "moccasin", "navajowhite", "navy",
    "oldlace", "olive", "olivedrab", "orange", "orangered",
    "orchid", "palegoldenrod", "palegreen", "paleturquoise",
    "palevioletred", "papayawhip", "peachpuff", "peru", "pink",
    "plum", "powderblue", "purple", "red", "rosybrown",
    "royalblue", "saddlebrown", "salmon", "sandybrown",
    "seagreen", "seashell", "sienna", "silver", "skyblue",
    "slateblue", "slategray", "slategrey", "snow", "springgreen",
    "steelblue", "tan", "teal", "thistle", "tomato", "turquoise",
    "violet", "wheat", "white", "whitesmoke", "yellow",
    "yellowgreen"]

    available_colorsequences = {"Plotly":px.colors.qualitative.Plotly, "G10":px.colors.qualitative.G10,
    "T10":px.colors.qualitative.T10, "Alphabet":px.colors.qualitative.Alphabet, "Dark24":px.colors.qualitative.Dark24, "Dark24_r":px.colors.qualitative.Dark24_r,
    "Light24":px.colors.qualitative.Light24, "Set1":px.colors.qualitative.Set1, "Pastel1":px.colors.qualitative.Pastel,
    "Dark2":px.colors.qualitative.Dark2, "Set2":px.colors.qualitative.Set2, "Pastel2":px.colors.qualitative.Pastel2,
    "Set3":px.colors.qualitative.Set3, "Antique":px.colors.qualitative.Antique,"Bold":px.colors.qualitative.Bold,
    "Pastel":px.colors.qualitative.Pastel, "Prism":px.colors.qualitative.Prism, "Safe":px.colors.qualitative.Safe,
    "Vivid":px.colors.qualitative.Vivid}

    st.session_state['available_colorsequences'] = available_colorsequences

    available_colorscales = px.colors.named_colorscales()

    available_taxonomic_levels_list= ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Hash']
    st.session_state['available_taxonomic_levels_list'] = available_taxonomic_levels_list

    return [directories_to_create, available_templates_list, available_clustering_units, plotly_colors, available_colorsequences, available_colorscales, available_taxonomic_levels_list]

directories_to_create, available_templates_list, available_clustering_units, plotly_colors, available_colorsequences, available_colorscales, available_taxonomic_levels_list = TTT_variables()

########################################################################################################################
# Collect path to default files
path_to_ttt = Path(__file__).resolve().parent
user_preferences_xlsx = path_to_ttt.joinpath('user_preferences.xlsx')

TTT_input_path_to_projects = ""
if user_preferences_xlsx.exists():
    user_data_df = pd.read_excel(user_preferences_xlsx).fillna('')
    TTT_input_path_to_projects = user_data_df[user_data_df['Variable'] == 'TTT_input_path_to_projects']['Value'].values.tolist()[0]
    TTT_height = user_data_df[user_data_df['Variable'] == 'TTT_height']['Value'].values.tolist()[0]
    TTT_width = user_data_df[user_data_df['Variable'] == 'TTT_width']['Value'].values.tolist()[0]
    TTT_showlegend = user_data_df[user_data_df['Variable'] == 'TTT_showlegend']['Value'].values.tolist()[0]
    TTT_template = user_data_df[user_data_df['Variable'] == 'TTT_template']['Value'].values.tolist()[0]
    TTT_fontsize = user_data_df[user_data_df['Variable'] == 'TTT_fontsize']['Value'].values.tolist()[0]
    TTT_hash = user_data_df[user_data_df['Variable'] == 'TTT_hash']['Value'].values.tolist()[0]
    TTT_scattersize = user_data_df[user_data_df['Variable'] == 'TTT_scattersize']['Value'].values.tolist()[0]
    TTT_linewidth = user_data_df[user_data_df['Variable'] == 'TTT_linewidth']['Value'].values.tolist()[0]
    TTT_jitter = user_data_df[user_data_df['Variable'] == 'TTT_jitter']['Value'].values.tolist()[0]
    TTT_color1 = user_data_df[user_data_df['Variable'] == 'TTT_color1']['Value'].values.tolist()[0]
    TTT_color2 = user_data_df[user_data_df['Variable'] == 'TTT_color2']['Value'].values.tolist()[0]
    TTT_color3 = user_data_df[user_data_df['Variable'] == 'TTT_color3']['Value'].values.tolist()[0]
    TTT_colorsequence = user_data_df[user_data_df['Variable'] == 'TTT_colorsequence']['Value'].values.tolist()[0]
    TTT_colorscale = user_data_df[user_data_df['Variable'] == 'TTT_colorscale']['Value'].values.tolist()[0]

########################################################################################################################

with st.sidebar:
    # Get user input
    #####################################
    st.write(""" ## Working directory """)

    TTT_input_path_to_projects = Path(st.text_input('📂 Enter Path to APSCALE Projects', value=TTT_input_path_to_projects))
    st.session_state['TTT_input_path_to_projects'] = TTT_input_path_to_projects
    path_to_projects = TTT_input_path_to_projects / 'TTT_projects'

    if not TTT_input_path_to_projects.exists() or TTT_input_path_to_projects == Path(''):
        st.error(f'Please enter PATH to a valid folder!')

    elif not path_to_projects.exists():
        st.warning(f'Could not find a "TTT_project" folder here: "{path_to_projects}"')
        if st.button('Create "TTT_projects" folder', width='stretch'):
            os.makedirs(path_to_projects, exist_ok=True)
            os.makedirs(path_to_projects / 'Default_project', exist_ok=True)
            os.makedirs(path_to_projects / 'Default_project' / 'Import', exist_ok=True)
            os.makedirs(path_to_projects / 'Default_project' / 'TaXon_tables', exist_ok=True)
            st.success(f'New project Default_project was created!')
            st.success(f'Created "TTT_project" folder at: "{path_to_projects}"')
            st.rerun()
    else:
        # Collect all projects
        projects_dict = {Path(i).name: Path(i) for i in sorted(glob.glob(str(path_to_projects / '*')))}

        if st.button("🔄 Refresh files and folders", width='stretch'):
            st.session_state.pop("df_taxon_table", None)
            st.session_state.pop("df_taxon_table_uniq", None)
            st.session_state.pop("df_metadata_table", None)
            st.session_state.pop("df_samples", None)
            st.session_state.pop("df_traits_table", None)
            st.session_state.pop("df_hashes", None)

        #####################################
        st.write(""" ## Project folders """)

        # Select project folder
        available_projects_list = list(projects_dict.keys())
        available_projects_list = st.selectbox(label='📂 Select a TTT project', options=available_projects_list, key='active_project', index=0)

        # Create a new project
        new_project_name = st.text_input(label='Enter name of new project', key='new_project_name')
        if st.button('⚙️ Create new TTT project', width='stretch'):
            if new_project_name != '':
                os.makedirs(path_to_projects / new_project_name, exist_ok=True)
                os.makedirs(path_to_projects / new_project_name / 'Import', exist_ok=True)
                os.makedirs(path_to_projects / new_project_name / 'TaXon_tables', exist_ok=True)
                st.success(f'New project {new_project_name} was created!')
                st.rerun()
            else:
                st.error('Please enter a valid project name!')

        # Select TTT project
        active_project_path = path_to_projects / st.session_state['active_project']
        st.session_state['active_project_path'] = active_project_path
        if st.button('🗂️ Open Active Project', width='stretch'):
            open_folder(active_project_path)

        # Select Table
        st.write(""" ## Taxon Tables """)

        taxon_tables_path = active_project_path / 'TaXon_tables' / '*.xlsx'
        available_taxontables_list = glob.glob(str(taxon_tables_path))

        if available_taxontables_list == []:
            st.info('Please first create a Taxon Table!')

        else:

            if st.button('🔄 Refresh Available Tables', width='stretch'):
                available_taxontables_list = glob.glob(str(taxon_tables_path))
            available_taxontables_dict = {Path(i).name:Path(i) for i in available_taxontables_list}
            st.selectbox(label='Select a TaXon table:', options=sorted(available_taxontables_dict.keys()), key='selected_table_xlsx', index=0)

            # Select Taxon Table
            active_taxontable_file = available_taxontables_dict[st.session_state['selected_table_xlsx']]

            # Update session state
            st.session_state['taxon_table_xlsx'] = Path(active_taxontable_file)
            if st.button('📥 Load Table', type='primary', width='stretch'):
                # Update table name for display
                st.session_state['table_display_name'] = Path(active_taxontable_file)
                # Taxon Table -> df_taxon_table
                df_taxon_table = pd.read_excel(active_taxontable_file, sheet_name='Taxon Table').fillna('')
                hash_loc = df_taxon_table.columns.tolist().index('Hash')
                species_loc = df_taxon_table.columns.tolist().index('Species')
                taxa_cols = df_taxon_table.columns.tolist()[hash_loc+1:species_loc+1]
                st.session_state['taxa_cols'] = taxa_cols
                taxon_string = (df_taxon_table[taxa_cols].astype(str).agg(";".join, axis=1))
                similarity_loc = df_taxon_table.columns.get_loc("Similarity") + 1
                if 'Taxon' not in df_taxon_table.columns.tolist():
                    df_taxon_table.insert(similarity_loc, "Taxon", taxon_string)
                st.session_state['df_taxon_table'] = df_taxon_table
                # Metadata Table -> df_metadata_table
                st.session_state['df_metadata_table'] = pd.read_excel(active_taxontable_file, sheet_name='Metadata Table').fillna('')
                # Samples -> df_samples
                seq_loc = st.session_state['df_taxon_table'].columns.tolist().index('Seq')
                taxontable_samples = st.session_state['df_taxon_table'].columns.tolist()[seq_loc+1:]
                metadata_samples = st.session_state['df_metadata_table']['Samples'].values.tolist()
                st.session_state['df_samples'] = metadata_samples
                species_loc = st.session_state['df_taxon_table'].columns.tolist().index('Species')
                # Traits -> df_traits_table
                st.session_state['df_traits_table'] = st.session_state['df_taxon_table'].columns.tolist()[species_loc+1:seq_loc]
                # Hashes -> df_hashes
                st.session_state['df_hashes'] = st.session_state['df_taxon_table']['Hash']
                # Taxon columns
                st.session_state['df_taxon_cols'] = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

                if check_reads() == False:
                    st.error('⚠️ Your read table seems to include non-numeric values! Please check your Taxon Table! ⚠️')
                    del st.session_state["df_taxon_table"]
                    del st.session_state["table_display_name"]

            if 'df_taxon_table' in st.session_state and 'table_display_name' in st.session_state:
                st.success(f'Active table: {st.session_state["table_display_name"].stem}')
                if st.button('📄 Open Active Table', width='stretch'):
                    open_file(st.session_state["table_display_name"])
            else:
                st.warning(f'No active table.')

            ################################################################################################################
            st.write(""" # Settings """)

            ################################################################################################################
            st.write(""" ### Presentation """)
            # PLOT HEIGHT
            st.session_state['TTT_height'] = st.sidebar.number_input('Height', 400, 4000, TTT_height)
            # PLOT WIDTH
            st.session_state['TTT_width'] = st.sidebar.number_input('Width', 400, 4000, TTT_width)
            # SHOW LEGEND
            st.session_state['TTT_showlegend'] = st.sidebar.selectbox('Show legend', [True, False], index=[True, False].index(TTT_showlegend))
            # TEMPLATE
            st.session_state['TTT_template'] = st.sidebar.selectbox('Layout', available_templates_list, index=available_templates_list.index(TTT_template))
            # FONT_SIZE
            st.session_state['TTT_fontsize'] = st.sidebar.slider('Font size', 6, 30, TTT_fontsize)
            # CLUSTERING UNIT
            st.session_state['TTT_hash'] = st.sidebar.selectbox('Clustering unit (hash description)', available_clustering_units, index=available_clustering_units.index(TTT_hash))
            # SCATTER SIZE
            st.session_state['TTT_scattersize'] = st.sidebar.number_input('Scatter size', 0, 40, TTT_scattersize)
            # LINE WIDTH
            st.session_state['TTT_linewidth'] = st.sidebar.number_input('Line width', 0, 40, TTT_linewidth)
            # JITTER
            st.session_state['TTT_jitter'] = st.sidebar.number_input('Jitter', 0.0, 1.0, 0.3)

            TTT_config = {
                "displaylogo": False,  # remove plotly logo
                "toImageButtonOptions": {
                    "format": "svg",  # png, svg, jpeg, webp
                    "filename": "TTT_plot",
                    "height": st.session_state['TTT_height'],
                    "width": st.session_state['TTT_width']
                },
                "modeBarButtonsToRemove": [
                    "lasso2d",
                    "select2d"
                ]}
            st.session_state['TTT_config'] = TTT_config

            ################################################################################################################
            st.sidebar.write(""" ### Colors """)
            # COLOR 1
            st.session_state['TTT_color1'] = st.sidebar.selectbox('Color #1', plotly_colors, index=plotly_colors.index(TTT_color1))
            # COLOR 2
            st.session_state['TTT_color2'] = st.sidebar.selectbox('Color #2', plotly_colors, index=plotly_colors.index(TTT_color2))
            # COLOR 3
            st.session_state['TTT_color3'] = st.sidebar.selectbox('Color #3', plotly_colors, index=plotly_colors.index(TTT_color3))

            # COLORSEQUENCE
            st.session_state['TTT_colorsequence'] = st.sidebar.selectbox('Color sequence', list(available_colorsequences.keys()), index=list(available_colorsequences.keys()).index(TTT_colorsequence))

            # COLORSCALE
            st.session_state['TTT_colorscale'] = st.sidebar.selectbox('Color scale', available_colorscales, index=available_colorscales.index(TTT_colorscale))

            settings_keys = [
                'TTT_input_path_to_projects',
                'TTT_height',
                'TTT_width',
                'TTT_showlegend',
                'TTT_template',
                'TTT_fontsize',
                'TTT_hash',
                'TTT_scattersize',
                'TTT_linewidth',
                'TTT_jitter',
                'TTT_color1',
                'TTT_color2',
                'TTT_color3',
                'TTT_colorsequence',
                'TTT_colorscale'
            ]

            settings_df = pd.DataFrame({
                "Variable": settings_keys,
                "Value": [st.session_state[key] for key in settings_keys]
            })

            if st.button("💾 Save Settings", width='stretch'):
                settings_df.to_excel(user_preferences_xlsx, index=False)
                st.success('Settings were saved!')

########################################################################################################################
if 'df_taxon_table' not in st.session_state or 'table_display_name' not in st.session_state:
    expanded=True
else:
    expanded=False

st.markdown("## 🐙 TTTutorial")
with st.expander(expanded=expanded, label='See more'):
    TTTutorial()

if path_to_projects.exists():

    st.markdown("## 🧬 Create Taxon Table")
    with st.expander(expanded=expanded, label='See more'):
        st.success(
            "This section allows you to generate a **Taxon Table** by combining read and taxonomy data from your project imports. "
            "Select the appropriate read and taxonomy tables, choose import formats, and provide a name for the new Taxon Table. "
            "The resulting table will be ready for downstream analyses such as diversity calculations, plotting, and ecological assessments."
        )

        # Collect import files
        import_folder = active_project_path / 'Import'
        import_folder_files = sorted([p for p in import_folder.glob('*') if p.suffix in ('.xlsx', '.snappy') and not p.name.startswith(('~$', '.'))])
        import_folder_files = {i.name:i for i in import_folder_files}
        st.session_state['import_folder_files'] = import_folder_files

        a1, a2 = st.columns(2)
        with a1:
            selected_read_table = st.selectbox(label='Select Read Table', options=import_folder_files.keys(), key='selected_read_table')
            st.selectbox(label='Choose Import Format', options=['APSCALE'], key='read_table_import_format')
        with a2:
            selected_taxonomy_table = st.selectbox(label='Select Taxonomy Table', options=import_folder_files.keys(), key='selected_taxonomy_table')
            st.selectbox(label='Choose Import Format', options=['APSCALE', 'BOLDigger'], key='taxonomy_table_import_format')

        st.text_input(label='Enter Name of Taxon Table', key='taxon_table_creation_name')

        if selected_read_table and selected_taxonomy_table:
            if st.button('Create Taxon Table'):
                create_taxon_table()

########################################################################################################################

if 'df_taxon_table' in st.session_state and 'table_display_name' in st.session_state:

    # general useful variables
    available_metadata = st.session_state['df_metadata_table'].columns.tolist()[1:]
    st.session_state['taxa_cols2'] = st.session_state['taxa_cols'] + ['Hash']

    ########################################################################################################################
    st.markdown("## 🖥️ Basic Stats")
    with st.expander(expanded=True, label='See more'):
        st.success(
            "This section provides an overview of key statistics describing your dataset. "
            "It summarizes sequencing depth, OTU richness, taxonomic resolution, and the most abundant taxa across samples. "
            "These metrics help assess data quality and give a first impression of biodiversity patterns in your dataset.")

        col1, col2 = st.columns(2)
        with col1:
            basic_stats_reads()
        with col2:
            basic_stats_OTUs()
        col1, col2 = st.columns(2)
        with col1:
            tax_res_plot()
        with col2:
            top_n_taxa_plot()
        try:
            st.write('Table 1: Dataset statistics')
            collect_sample_stats()
        except:
            st.warning('Unable to calculate dataset statistics!')

    ########################################################################################################################
    st.markdown("## 🛠️ Table Processing")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "This section allows you to process your taxon table before downstream analyses. "
            "You can merge replicates, subtract negative controls, filter by reads, samples, taxa, or traits, "
            "and apply normalization to make your data comparable across samples. ")

        st.write('### Replicate Merging')
         # --- Overview Metrics ---
        col1, col2, col3, = st.columns(3)
        with col1:
            suffix_detection = st.selectbox(label='Replicate Suffixes', options=['Auto-detect', 'Custom'], key='suffix_detection')
            if suffix_detection == 'Auto-detect':
                # --- Derive sample structure ---
                df_samples = st.session_state['df_samples']
                hash_label = st.session_state['TTT_hash']
                suffixes = sorted({i.split('_')[-1] for i in df_samples})
                unique_samples = sorted({'_'.join(i.split('_')[:-1]) for i in df_samples})
                st.session_state['suffixes'] = suffixes
                st.session_state['unique_samples'] = unique_samples
                expected_samples = [f"{s}_{suf}" for suf in suffixes for s in unique_samples]
                missing_samples = [sample for sample in expected_samples if sample not in df_samples]
                n_suffixes = len(suffixes)
                cutoff_options = list(range(1, n_suffixes + 1))[::-1]
                # --- Compact UI ---
                st.write(
                    f"""
                    **Detected suffixes:** {', '.join(suffixes)}  
                    **Unique samples:** {len(unique_samples)}  
                    **Expected replicates:** {len(expected_samples)}  
                    **Present replicates:** {len(df_samples)}  
                    **Missing samples:** {len(missing_samples)}  
                    """)
            else:
                # --- Derive sample structure ---
                df_samples = st.session_state['df_samples']
                hash_label = st.session_state['TTT_hash']
                suffixes_input = st.text_input(label='Enter suffixes:', value='A,B', key='suffixes_input')
                suffixes = suffixes_input.split(',')
                st.session_state['suffixes'] = suffixes
                unique_samples = sorted({'_'.join(i.split('_')[:-1]) for i in df_samples})
                st.session_state['unique_samples'] = unique_samples
                expected_samples = [f"{s}_{suf}" for suf in suffixes for s in unique_samples]
                missing_samples = [sample for sample in expected_samples if sample not in df_samples]
                n_suffixes = len(suffixes)
                cutoff_options = list(range(1, n_suffixes + 1))[::-1]
                # --- Compact UI ---
                st.write(
                    f"""
                    **Detected suffixes:** {', '.join(suffixes)}  
                    **Unique samples:** {len(unique_samples)}  
                    **Expected replicates:** {len(expected_samples)}  
                    **Present replicates:** {len(df_samples)}  
                    **Missing samples:** {len(missing_samples)}  
                    """)
        with col2:
            pass
        with col2:
            cutoff_value = st.selectbox(label=f"{hash_label} must be present in at least:",options=cutoff_options, format_func=lambda x: f"{x} replicate(s)", key='cutoff_value')
            if cutoff_value == 1:
                st.info(f"All {hash_label} will be kept after replicate merging.")
            else:
                st.warning(
                    f"Only {hash_label} present in "f"{cutoff_value}/{n_suffixes} replicates "f"will be kept after replicate merging.")
        with col3:
            missing_replicates_handling = st.selectbox(label=f"Handle samples which miss replicates:",options=['Keep', 'Remove'], key='missing_replicates_handling')
        # --- Action Button ---
        if st.button("🔄 Merge Replicates", width='stretch'):
            replicate_merging()
        st.divider()


        st.write('### Negative Control Subtraction')
        # --- Initialize structure ---
        col1, col2 = st.columns([4, 6])
        with col1:
            NC_identifier = st.text_input(label='Negative Control identifier string:', value='NC_', key='NC_identifier')
        with col2:
            NC_handling = st.selectbox(label='Negative Control handling behaviour:', options=['Subtract Sum', 'Subtract Max', 'Subtract Mean', f'Remove {st.session_state["TTT_hash"]}'], key='NC_handling')
        NC_default = [i for i in st.session_state['df_samples'] if NC_identifier in i]
        df_NC_samples = st.multiselect(label='Select Negative Controls:', options=st.session_state['df_samples'], default=NC_default)
        st.session_state['df_NC_samples'] = df_NC_samples
        st.info(f'Found {len(df_NC_samples)} Negative Controls in {len(st.session_state["df_samples"])} samples!')
        # --- Action Button ---
        if st.button("🔄 Subtract Negative Controls", width='stretch'):
            NC_subtraction()
        st.divider()


        st.write('### Read-based Filter')
        # --- Derive structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            read_filter = st.selectbox(label=f'Filter {st.session_state["TTT_hash"]} based on:', options=['Absolute Reads', 'Relative Reads'], key='read_filter')
        with col2:
            read_filter_mode = st.selectbox(label='Read Filter Mode:', options=['Per Sample (column)', f'Per {st.session_state["TTT_hash"][:-1]} (row)'], key='read_filter_mode')
        with col3:
            if read_filter == 'Absolute Reads':
                read_filter_value = st.number_input(label='Absolute Read Threshold:', value=10, key='read_filter_value')
            if read_filter == 'Relative Reads':
                read_filter_value = st.number_input(label='Relative Read Threshold (%):', value=1.00, key='read_filter_value', max_value=100.00)
        # --- Action Button ---
        if st.button("🔄 Apply read-based filter", width='stretch'):
            read_based_filter()
        st.divider()


        st.write('### Read-based Normalisation')
        # --- Derive structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            normalisation_mode = st.selectbox(label='Choose Normalisation Mode:', options=['Custom Read Threshold', 'Auto-detect Minimum Threshold'], key='normalisation_mode')
        with col2:
            if normalisation_mode == 'Auto-detect Minimum Threshold':
                st.session_state['sub_sample_size'] = st.session_state['df_taxon_table'][st.session_state['df_samples']].sum().min()
            elif normalisation_mode == 'Custom Read Threshold':
                st.session_state['sub_sample_size'] = st.number_input(label='Sub-sample Read Threshold:', value=250000, min_value=1)
        with col3:
            mean_reads = round(st.session_state['df_taxon_table'][st.session_state['df_samples']].sum().mean(),2)
            samples_below = [i for i in st.session_state['df_taxon_table'][st.session_state['df_samples']].sum().values.tolist() if i >= st.session_state['sub_sample_size']]
            st.info(f'Mean: {mean_reads:,} reads per sample.')
            st.info(f"{len(samples_below)}/{len(st.session_state['df_samples'])} samples fall within the normalisation threshold.")
        st.info(f"Normalsation value: {st.session_state['sub_sample_size']:,} reads.")
        # --- Action Button ---
        if st.button("🔄 Apply read-based normalisation", width='stretch'):
            read_based_normalisation()
        st.divider()


        st.write('### Taxonomic Filter')
        col1, col2 = st.columns([2,6])
        with col1:
            taxon_select = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols'], key='taxon_select')
            taxon_filter_action = st.selectbox(label='Select an Action:', options=['Keep', 'Remove'], key='taxon_filter_action')
            taxon_filter_suffix = st.text_input(label='Enter Table Suffix:', value=f'{taxon_select.lower()}')
            st.session_state['taxon_filter_suffix'] = taxon_filter_suffix
        with col2:
            available_taxa = sorted([i for i in st.session_state['df_taxon_table'][taxon_select].unique() if i != ''])
            selected_taxa = st.multiselect(label='Select Taxa:', options=available_taxa, default=available_taxa, key='selected_taxa')
        st.info(f"Selected {len(selected_taxa)} of {len(available_taxa)} taxa.")
        # --- Action Button ---
        if st.button("🔄 Apply taxonomy filter", width='stretch'):
            taxonomic_filter()
        st.divider()


        st.write('### Sample Filter')
        col1, col2 = st.columns([2,6])
        with col1:
            sample_selector_str = st.text_input(label='Enter STRING to Filter Samples:', value='', key='sample_selector_str')
            sample_filter_action = st.selectbox(label='Select an Action:', options=['Keep', 'Remove'], key='sample_filter_action')
            sample_filter_suffix = st.text_input(label='Enter Table Suffix:', value='sample_filter', key='sample_filter_suffix')
        with col2:
            available_samples = sorted(st.session_state['df_samples'])
            selected_samples = [i for i in available_samples if sample_selector_str in i]
            samples_to_filter = st.multiselect(label='Select Samples:', options=available_samples, default=selected_samples)
            st.session_state['samples_to_filter'] = samples_to_filter
        st.info(f"Selected {len(selected_samples)} of {len(available_samples)} samples.")
        # --- Action Button ---
        if st.button("🔄 Apply sample filter", width='stretch'):
            sample_filter()
        st.divider()


        st.write('### Trait Filter')
        col1, col2 = st.columns([2,6])
        with col1:
            available_traits = sorted(st.session_state['df_traits_table'])
            selected_trait = st.selectbox(label='Select Trait:', options=available_traits)
            st.session_state['selected_trait'] = selected_trait
            trait_filter_action = st.selectbox(label='Select an Action:', options=['Keep', 'Remove'], key='trait_filter_action')
            trait_filter_suffix = st.text_input(label='Enter Table Suffix:', value=selected_trait)
            st.session_state['trait_filter_suffix'] = trait_filter_suffix
        with col2:
            available_trait_values = sorted([str(i) for i in set(st.session_state['df_taxon_table'][selected_trait])])
            series = pd.Series(available_trait_values)
            numeric = pd.to_numeric(series, errors="coerce")
            if numeric.notna().all():
                # All values are numeric
                available_trait_values = numeric.tolist()
                trait_type = "numeric"
            else:
                # Not all numeric → convert everything to string
                available_trait_values = series.astype(str).tolist()
                trait_type = "string"
            if trait_type == "numeric":
                numeric_values = [float(i) for i in available_trait_values]
                selected_trait_values = st.number_input(label='Select Threshold:', value=max(numeric_values), max_value=max(numeric_values), min_value=min(numeric_values))
                st.session_state['selected_trait_values'] = selected_trait_values
                st.session_state['trait_type'] = 'numeric'
                if trait_filter_action == 'Keep':
                    st.info(f"All {st.session_state['TTT_hash']} with '{selected_trait} >= {selected_trait_values}' will be kept!")
                else:
                    st.info(f"All {st.session_state['TTT_hash']} with '{selected_trait} <= {selected_trait_values}' will be kept!")
            else:
                selected_trait_values = st.multiselect(label='Select Trait Values:', options=available_trait_values, default=available_trait_values)
                st.session_state['trait_type'] = 'string'
                st.info(f"Selected {len(selected_trait_values)} of {len(available_trait_values)} trait values.")
        # --- Action Button ---
        if st.button("🔄 Apply trait filter", width='stretch'):
            trait_filter()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🔁 Table Conversion")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "This section provides tools to convert and enhance your taxon tables. "
            "You can simplify redundant entries, combine multiple tables, import traits, sort or rename samples, "
            "and export your data to FASTA format. Use the buttons below to apply each step interactively.")

        st.write('### Simplify Taxon Table')
        # --- Initialize structure ---
        st.info(
            f"{st.session_state['TTT_hash']} with identical taxonomic assignments "
            f"will be merged into a single entry. The representative "
            f"{st.session_state['TTT_hash']} will be selected based on the highest "
            f"read abundance and sequence similarity."
        )
        # --- Action Button ---
        if st.button("🔄 Simplify Taxon Table", width='stretch'):
            simplify_table()
        st.divider()

        st.write('### Combine Taxon Tables')
        # --- Initialize structure ---
        col1, col2 = st.columns([4, 6])
        with col1:
            df2_taxon_table_file = st.file_uploader(label='Drag and drop your trait table here:', key='df2_taxon_table_file')
        with col2:
            table_merge_suffix = st.text_input(label='Please Enter a Table Suffix:', value='combined', key='table_merge_suffix')
            if df2_taxon_table_file == None:
                st.write('')
                st.info('Please provide a Taxon Table to continue!')
        # --- Action Button ---
        if st.button("🔄 Merge Selected Taxon Tables", width='stretch'):
            merge_taxon_tables()
        st.divider()

        st.write('### Add Traits From File')
        # --- Initialize structure ---
        col1, col2 = st.columns([4, 6])
        with col1:
            trait_import_table = st.file_uploader(label='Drag and drop your trait table here:', key='trait_import_table')
        with col2:
            if trait_import_table == None:
                st.write('')
                st.info('Please provide an Excel table (.xlsx) to continue!')
            elif not str(trait_import_table.name).endswith('.xlsx'):
                st.write('')
                st.error('Please upload an Excel file (.xlsx)')
            else:
                trait_import_df = pd.read_excel(trait_import_table)
                st.session_state['trait_import_df'] = trait_import_df
                st.write('')
                st.info(f'Found Dataframe with {len(trait_import_df)} taxa and {len(trait_import_df.columns)-1} traits!')
                st.session_state['trait_import_taxon'] = trait_import_df.columns.tolist()[0]
                st.info(f'Selected taxon: {st.session_state["trait_import_taxon"]}')
        # --- Action Button ---
        if st.button("🔄 Add Traits From File", width='stretch'):
            if trait_import_table != None:
                add_traits_from_file()
        st.divider()


        st.write('### Sort Samples')
        # --- Initialize structure ---
        st.info('Sort your Metadata Table as desired. The samples of the Taxon Table will then be sorted accordingly!')
        # --- Action Button ---
        if st.button("🔄 Sort Samples", width='stretch'):
            sort_samples()
        st.divider()


        st.write('### Rename Samples')
        # --- Initialize structure ---
        col1, col2 = st.columns([4, 6])
        with col1:
            sample_rename_metadata = st.selectbox(label='Select Metadata Category', options=available_metadata, key='sample_rename_metadata')
        with col2:
            st.write('')
            st.info(f'Samples will be renamed according to metadata category "{sample_rename_metadata}". Old samples names are stored in the Metadata Table.')
        # --- Action Button ---
        if st.button("🔄 Rename Samples", width='stretch'):
            sample_rename()
        st.divider()


        st.write('### Export Table To Fasta')
        # --- Initialize structure ---
        st.info(f"{st.session_state['TTT_hash']} will be exported to .FASTA format.")
        # --- Action Button ---
        if st.button("🔄 Export Table To Fasta", width='stretch'):
            export_fasta()
        st.divider()

    ########################################################################################################################
    st.markdown("## 📈 Metabarcoding Basics")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "This section allows you to explore the distribution and diversity of reads across your samples. "
            "You can generate read distribution diagrams, perform read-based rarefaction to assess sampling depth, "
            "and check read-to-OTU/sequence autocorrelation per sample."
        )

        st.write(f'### Read/{st.session_state["TTT_hash"][:-1]} Distribution Plot')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            readdist_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols'], index=2, key='readdist_taxon')
            readdist_level = st.selectbox(label='Select Sample Level:', options=['Samples'] + available_metadata, key='readdist_level')
        with col2:
            readdist_mode = st.selectbox(label='Select Display Mode:', options=['Relative Reads', 'Absolute Reads', st.session_state["TTT_hash"]], key='readdist_mode')
            readdist_nan = st.selectbox(label='Handle Missing Taxonomy:', options=['Exclude', 'Include'], key='readdist_nan')
        with col3:
            readdist_color = st.selectbox(label='Select Color Mode:', options=['Color Sequence', 'Color Scale'], key='readdist_color')
        # --- Action Button ---
        if st.button(f"🔄 Generate Read Distribution Diagram", width='stretch'):
            readdist_diagram()
        st.divider()

        st.write('### Circle Diagram')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)

        with col1:
            pycircos_taxon0 = st.selectbox(label='Select Outer Groups:', options=st.session_state['taxa_cols2'][:-3], index=0, key='pycircos_taxon0')
            available_taxon1 = st.session_state['taxa_cols2'][st.session_state['taxa_cols2'].index(pycircos_taxon0)+1:]
            pycircos_taxon1 = st.selectbox(label='Select Inner Groups:', options=available_taxon1, index=0, key='pycircos_taxon1')
            available_taxon2 = st.session_state['taxa_cols2'][st.session_state['taxa_cols2'].index(pycircos_taxon1)+1:]
            pycircos_taxon2 = st.selectbox(label='Alpha Diversity Measure:', options=available_taxon2, index=available_taxon2.index('Species'), key='pycircos_taxon2')
        with col2:
            pycircos_groups = [i for i in st.session_state['df_taxon_table'][pycircos_taxon0].unique().tolist() if i != '']
            st.session_state['pycircos_groups'] = pycircos_groups
            n_groups = len(pycircos_groups)
            pycircos_cmaps_dict = {'Green': 'Greens', 'Orange': 'Oranges', 'Yellow': 'Wistia', 'Red': 'Reds',
                                   'Blue': 'Blues', 'Lightgrey': 'Greys'}
            if n_groups <= 6:
                st.session_state['pycircos_colors'] = 'custom'
                st.session_state['pycircos_cmaps_dict'] = pycircos_cmaps_dict
                cols = st.columns(2)
                for i, group in enumerate(pycircos_groups):
                    with cols[i % 2]:
                        st.selectbox(
                            label=group,
                            options=list(pycircos_cmaps_dict.keys()),
                            index=i,
                            key=f'pycircos_colors_{group}')
            else:
                pycircos_colors = st.selectbox(label='Select Color Mode:', options=list(pycircos_cmaps_dict.keys()), key='pycircos_colors')
        with col3:
            pycircos_nan = st.selectbox(label='Handle Missing Taxonomy:', options=['Exclude', 'Include'], key='pycircos_nan')
            pycircos_inner_labels = st.selectbox(label='Inner Labels:', options=['Show', 'Hide'], key='pycircos_inner_labels')
        # --- Action Button ---
        if st.button(f"🔄 Generate Pycircos Diagram", width='stretch'):
            pycircos_plot()
        st.divider()

        st.write(f'### Read-based Rarefaction')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            read_rarefaction_reps = st.number_input(label='Number of Repetitions:', min_value=10, max_value=5000, value=10, key='read_rarefaction_reps')
            read_rarefaction_splits = st.number_input(label='Number of Sampling Points (X-Axis):', min_value=10, max_value=100, value=10, key='read_rarefaction_splits')
        with col2:
            read_rarefaction_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols2'], index=6, key='read_rarefaction_taxon')
            read_rarefaction_display = st.selectbox(label='Select Y-Axis Scaling:', options=['Absolute', 'Relative'], key='read_rarefaction_display')
        with col3:
            read_rarefaction_color = st.selectbox(label='Select Color Mode:', options=['Color Scale', 'Color Sequence', 'Single Color'], index=2, key='read_rarefaction_color')
        st.info('Please note that read-based rarefaction analyses can take very long, especially with high sequencing depths!')
        # --- Action Button ---
        if st.button(f"🔄 Calculate Rarefaction Curves", width='stretch'):
            read_rarefaction()
        st.divider()

        st.write(f'### Read/{TTT_hash[:-1]} Auto-correlation')
        # --- Initialize structure ---
        st.info(f'Correlate the number of Reads and {TTT_hash} per sample.')
        # --- Action Button ---
        if st.button(f"🔄 Calculate Read/{TTT_hash[:-1]} Auto-correlation", width='stretch'):
            read_hash_autocorrelation()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🌱 Alpha Diversity")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "Explore the alpha diversity of your samples in this section. "
            "You can generate alpha diversity diagrams using various metrics and visualization types, "
            "and visualize hierarchical relationships with PyCircos circle diagrams."
        )

        st.write('### Alpha Diversity Diagram')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            alphadiv_taxon = st.selectbox(label='Alpha Diversity Measure:', options=st.session_state['taxa_cols2'], index=6, key='alphadiv_taxon')
            alphadiv_measure = st.selectbox(label='Alpha Diversity Calculation:', options=['Richness', 'Heip Evenness', 'Pielou Evenness', 'Shannon Diversity'], key='alphadiv_measure', help=help_alphadiv)
            if alphadiv_measure != 'Richness' and alphadiv_taxon == 'Hash':
                alphadiv_mode = st.selectbox(label='Select Quantification:', options=['Rel. Reads'], key='alphadiv_mode')
            elif alphadiv_measure != 'Richness' and alphadiv_taxon != 'Hash':
                alphadiv_mode = st.selectbox(label='Select Quantification:', options=[f'{TTT_hash[:-1]} Richness', 'Rel. Reads'], key='alphadiv_mode')
            else:
                alphadiv_mode = st.selectbox(label='Select Quantification:', options=['None'], key='alphadiv_mode')
        with col2:
            alphadiv_layout = st.selectbox(label='Select Display Mode:', options=['Bar', 'Scatter', 'Box', 'Violin'], key='alphadiv_layout')
            if alphadiv_layout == 'Bar' or alphadiv_layout == 'Scatter':
                alphadiv_level = st.selectbox(label='Select Sample Level:', options=['Samples'], key='alphadiv_level')
            else:
                alphadiv_level = st.selectbox(label='Select Sample Level:', options=available_metadata, key='alphadiv_level')
            if alphadiv_measure != 'Richness' and alphadiv_measure != 'Shannon Diversity':
                alphadiv_interpration = st.selectbox(label='Add Interpretation Lines:', options=['Yes', 'No'], key='alphadiv_interpration')
            else:
                st.session_state['alphadiv_interpration'] = 'No'
        with col3:
            alphadiv_color = st.selectbox(label='Select Color Mode:', options=['Color Scale', 'Color Sequence', 'Single Color'], key='alphadiv_color')
            if alphadiv_layout == 'Box' or alphadiv_layout == 'Violin':
                alphadiv_boxpoints = st.selectbox(label='Boxpoints Display:', options=['All', 'False'], key='alphadiv_boxpoints')
            else:
                st.session_state['alphadiv_boxpoints'] = 'NAN'
        # --- Action Button ---
        if st.button(f"🔄 Generate {alphadiv_taxon} {alphadiv_measure} Diagram", width='stretch'):
            alphadiv_diagram()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🌎 Beta Diversity")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "This section allows you to explore beta diversity across your samples. "
            "You can generate distance/dissimilarity matrices, perform PCoA and NMDS analyses, "
            "and visualize multivariate relationships between samples and groups."
        )

        st.write('### Beta Diversity Matrix')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            betadiv_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols2'], index=6, key='betadiv_taxon')
            betadiv_level = st.selectbox(label='Select Sample Level:', options=['Samples'] + available_metadata, key='betadiv_level')
        with col2:
            betadiv_calculation = st.selectbox(label='Select Diversity Calculation:', options=['Jaccard (qualitative)', 'Bray-Curtis (quantitative)'], key='betadiv_mode')
            betadiv_measure = st.selectbox(label='Select Diversity Measure:', options=['Similarity', 'Dissimilarity'], key='betadiv_measure')
        with col3:
            betadiv_color = st.selectbox(label='Matrix Values:', options=['Show', 'Hide'], key='betadiv_labels')
            richness_nan = st.selectbox(label='Handle Missing Taxonomy:', options=['Exclude', 'Include'], key='betadiv_nan')
        # --- Action Button ---
        if st.button(f"🔄 Generate {betadiv_calculation} {betadiv_measure} Diagram", width='stretch'):
            beta_diversity_matrix()
        st.divider()

        st.write('### Principle Coordinate (PCoA) Analysis')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            pcoa_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols2'], index=6, key='pcoa_taxon')
            pcoa_mode = st.selectbox(label='Select Diversity Calculation:', options=['Jaccard (qualitative)', 'Bray-Curtis (quantitative)'], key='pcoa_mode')
        with col2:
            pcoa_level = st.selectbox(label='Select Group Level:', options=['Samples'] + available_metadata, key='pcoa_level')
            pcoa_showlabels = st.selectbox(label='Scatter Labels:', options=showlabels_options, key='pcoa_showlabels')
        with col3:
            pcoa_color = st.selectbox(label='Select Color Mode:', options=['Color scale', 'Color sequence', 'Single Color'], key='pcoa_color')
            pcoa_groups = st.selectbox(label='Select Grouping Display:', options=['None', 'Outlines', 'Continous'], key='pcoa_groups')
        # --- Action Button 1 ---
        if st.button(f"🔄 Calculate {pcoa_mode} Dissimilarity PCoA", width='stretch'):
            pcoa_df, explained_variance_df, pcoa_distance_matrix_df = pcoa_calculation()
            st.session_state['pcoa_df'] = pcoa_df
            st.session_state['pcoa_distance_matrix_df'] = pcoa_distance_matrix_df
            st.session_state['pcoa_explained_variance_df'] = explained_variance_df
            st.session_state['pcoa_explained_variance_dict'] = {f'{i} ({round(j[0],2)}%)':i for i,j in zip(explained_variance_df.index.tolist(), explained_variance_df.values.tolist())}
        # --- Action Button 2 ---
        if 'pcoa_explained_variance_dict' in st.session_state and 'pcoa_df' in st.session_state:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.selectbox(label='Select X-Axis:', options=st.session_state['pcoa_explained_variance_dict'].keys(), index=0, key='pcoa_x')
            with col2:
                st.selectbox(label='Select Y-Axis:', options=st.session_state['pcoa_explained_variance_dict'].keys(), index=1, key='pcoa_y')
            with col3:
                st.selectbox(label='Select Z-Axis:', options=['None'] + list(st.session_state['pcoa_explained_variance_dict'].keys()), index=0, key='pcoa_z')
            if st.button(f"🔄 Plot PCoA", width='stretch'):
                pcoa_plot()
        else:
            st.warning('Please calculate the PCoA first!')
        st.divider()

        st.write('### Non-metric Multidimensional Scaling (NMDS) Analysis')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            nmds_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols2'], index=6, key='nmds_taxon')
            nmds_mode = st.selectbox(label='Select Diversity Calculation:', options=['Jaccard (qualitative)', 'Bray-Curtis (quantitative)'], key='nmds_mode')
        with col2:
            nmds_level = st.selectbox(label='Select Group Level:', options=['Samples'] + available_metadata, key='nmds_level')
            nmds_showlabels = st.selectbox(label='Scatter Labels:', options=showlabels_options, key='nmds_showlabels')
        with col3:
            nmds_color = st.selectbox(label='Select Color Mode:', options=['Color scale', 'Color sequence', 'Single Color'], key='nmds_color')
            nmds_groups = st.selectbox(label='Select Grouping Display:', options=['None', 'Outlines', 'Continous'], key='nmds_groups')
            nmds_axes = st.selectbox(label='Select NMDS Dimensions:', options=[2,3], index=0, key='nmds_axes')
        # --- Action Button 1 ---
        if st.button(f"🔄 Calculate {nmds_mode} Dissimilarity NMDS", width='stretch'):
            nmds_df, nmds_stress = nmds_calculation()
            st.session_state['nmds_df'] = nmds_df
            st.session_state['nmds_stress'] = nmds_stress
        # --- Action Button 2 ---
        if 'nmds_stress' in st.session_state:
            st.info(f'NMDS stress: {st.session_state["nmds_stress"]:.2}')
            if st.button(f"🔄 Plot NMDS", width='stretch'):
                nmds_plot()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🧪 Sample Comparison")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "Use this section to compare your samples based on metadata categories. "
            "You can generate Venn diagrams to visualize shared and unique taxa between 2–3 groups. "
            "For more than 3 categories, consider using an UpSet chart instead."
        )

        st.write('### Venn Diagram')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            venn_metadata = st.selectbox(label='Select Metadata Category:', options=available_metadata, key='venn_metadata')
            venn_categories = [str(i) for i in st.session_state['df_metadata_table'][venn_metadata].unique() if i != '']
            st.info(f"Groups: {', '.join(list(venn_categories)[:10])}")
        with col2:
            venn_taxon_level = st.selectbox(label='Select taxonomic level:', options=st.session_state['taxa_cols2'], key='venn_taxon_level')
        with col3:
            venn_weighted = st.selectbox(label='Weighted Circle Sizes:', options=['Weighted', 'Non-weighted'], key='venn_weighted')
            venn_display = st.selectbox(label='Display Values:', options=['Values', 'Rel. Values', 'Both'], key='venn_display')
        # --- Action Button ---
        if len(venn_categories) >= 4:
            st.error('Please choose a Metadata with 2-3 Categories! Use an UpSet Chart instead!')
        elif st.button(f"🔄 Generate Venn{len(venn_categories)} Diagram", width='stretch'):
            venn_diagram()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🔍 Population Dynamics")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "Analyze population-level variation in your metabarcoding dataset. "
            "Haplotype distribution plots allow you to examine how genetic variants within a taxon "
            "are distributed across samples, providing insights into population structure and diversity."
        )

        st.write('### Haplotype Distribution Plot')
        col1, col2, col3 = st.columns(3)
        with col1:
            haplotype_level = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols'], index=6, key='haplotype_level')
        with col2:
            available_taxa = sorted([i for i in st.session_state['df_taxon_table'][haplotype_level].unique() if i != ''])
            haplotype_taxa = st.multiselect(label='Select Taxa:', options=available_taxa, default=available_taxa, key='haplotype_taxa')
        with col3:
            haplotype_mode = st.selectbox(label='Select Display Mode:', options=['Bar', 'Scatter'], key='haplotype_mode')

        if st.button(f"🔄 Generate Haplotype Distribution Plot", width='stretch'):
            haplotype_distribution()
        st.divider()

    ########################################################################################################################
    st.markdown("## 📍 Biogeography")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "Explore the biogeographic patterns of your metabarcoding data by visualizing sampling sites and taxon distributions on interactive maps. "
            "The GBIF Occurrence Validation compares detected taxa with known occurrence records within a defined radius around each sampling site. "
            "You can also generate distribution maps to examine how reads or relative abundances of selected taxa vary across sampling locations."
        )

        cols = st.session_state['df_metadata_table'].columns.tolist()
        if "Latitude" in cols and "Longitude" in cols:
            st.write('### GBIF Occurrence Validation')
            col1, col2 ,col3 = st.columns(3)
            with col1:
                available_map_styles = ["outdoors", "open-street-map", "carto-positron", "carto-darkmatter"]
                map_style = st.selectbox(label='Select Map Style:', options=available_map_styles, index=2, key='map_style')
                map_save = st.selectbox(label='Save Map as PDF and HTML:', options=['Yes', 'No'], index=1, key='map_save')
            with col2:
                map_color = st.selectbox(label='Select Color Mode:', options=['Color Scale', 'Color Sequence', 'Single Color'], index=2, key='map_color')
                map_groups = st.selectbox(label='Select Color Group:', options=['Samples'] + available_metadata, key='map_groups')
            with col3:
                sample_radius = st.number_input(label='Select Sample Radius:', value=20, min_value=0, max_value=1000, key='sample_radius')
                map_zoom = st.number_input(label='Select Map Zoom:', value=6, min_value=1, max_value=15, key='map_zoom')
            if st.button(f"🔄 Run GBIF Occurrence Validation", width='stretch'):
                gbif_occurrence_validation()
            sample_map()
            st.divider()

            st.write('### Taxon Distribution Maps')
            col1, col2 ,col3 = st.columns(3)
            with col1:
                available_map_styles = ["outdoors", "open-street-map", "carto-positron", "carto-darkmatter"]
                distmap_style = st.selectbox(label='Select Map Style:', options=available_map_styles, index=2, key='distmap_style')
                distmap_save = st.selectbox(label='Save Map as PDF and HTML:', options=['Yes', 'No'], index=1, key='distmap_save')
            with col2:
                distmap_level = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols'], index=6, key='distmap_level')
                available_taxa = sorted([i for i in st.session_state['df_taxon_table'][distmap_level].unique() if i != ''])
                distmap_taxon = st.selectbox(label='Select Taxon:', options=available_taxa, key='distmap_taxon')
            with col3:
                distmap_measure = st.selectbox(label='Select Measure:', options=['Rel. Reads', st.session_state['TTT_hash']], key='distmap_measure')
                distmap_zoom = st.number_input(label='Select Map Zoom:', value=6, min_value=1, max_value=15, key='distmap_zoom')
            taxon_distribution_map()
            st.divider()

        else:
            st.info('Please add the columns "Latitude" and "Longitude" to your Metadata Table!')

    ########################################################################################################################
    st.markdown("## 📆 Time Series")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "Visualize how taxon richness changes over time or across sample groups. "
            "This time series analysis allows you to assess trends, detect seasonal patterns, or monitor experimental changes. "
            "Confidence intervals and bootstrap repetitions can be adjusted to quantify uncertainty in richness estimates."
        )

        st.write('### Richness Over Time')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            ts_level = st.selectbox(label='Select Group Level:', options=['Samples'] + available_metadata, key='ts_level')
            ts_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols2'], index=6, key='ts_taxon')
        with col2:
            ts_ci = st.number_input('Select Confidence Intervall:', 0, 100, 98, key='ts_ci')
        with col3:
            ts_bootstraps = st.number_input('Select Bootstrap Repetitions:', 0, 5000, 100, key='ts_bootstraps')

        st.info(f'Samples will be sorted on the X-axis according to the metadata table "Samples" column order: {",".join(df_samples[:3])} ... ')
        # --- Action Button 1 ---
        if st.button(f"🔄 Calculate Richness Over Time", width='stretch'):
            time_series_richness()

    ########################################################################################################################
    st.markdown("## 🌿 GBIF Modules")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "The GBIF modules allow you to integrate your metabarcoding dataset with GBIF species data and tools. "
            "You can download occurrence records, enrich your taxon table, and prepare data for submission to GBIF. "
            "Ensure that your dataset meets the listed requirements to guarantee accurate processing and compatibility with GBIF workflows."
        )

        st.write('### GBIF Data')
        st.info("Download GBIF species data and enrich the taxon table.")
        if st.button("🔄 Add GBIF Data", width='stretch'):
            gbif_accession()
        st.divider()

        st.write('### GBIF Metabarcoding Data Toolkit')
        col1, col2, col3 = st.columns(3)
        with col1:
            gbif_target_gene = st.text_input('Target Gene:', value='', key='gbif_target_gene')
            gbif_date = st.selectbox(label='Select Collection Date Column from Metadata:', options=['None'] + available_metadata, key='gbif_date')
            st.write('')
            gbif_dataset_requirements = st.checkbox(label='**I have read the dataset requirements**.', value=False,
                                                    key='gbif_dataset_requirements')
        with col2:
            gbif_fwd_primer_name = st.text_input('Forward Primer Name:', value='', key='gbif_fwd_primer_name')
            gbif_fwd_primer_seq = st.text_input('Forward Primer Seq:', value='', key='gbif_fwd_primer_seq')
            gbif_lat = st.selectbox(label='Select Latitiude Column from Metadata:', options=['None'] + available_metadata, key='gbif_lat')
        with col3:
            gbif_rvs_primer = st.text_input('Forward Primer Name:', value='', key='gbif_rvs_primer_name')
            gbif_rvs_primer_seq = st.text_input('Forward Primer Seq:', value='', key='gbif_rvs_primer_seq')
            gbif_lon = st.selectbox(label='Select Longitude Column from Metadata:', options=['None'] + available_metadata, key='gbif_lon')

        if gbif_dataset_requirements == False:
            st.markdown("""
            #### 📋 Dataset Requirements
            Please ensure that the following conditions apply to your dataset before proceeding:
            - ✅ **DNA metabarcoding data**  
              The dataset originates from a DNA metabarcoding study.
            - ✅ **Data sharing permission**  
              You have permission to share and process the dataset.
            - ✅ **OTU/ASV table available**  
              The dataset contains an OTU/ASV table with associated taxonomy and sample metadata.
            - ✅ **Sequences available**  
              OTU/ASV sequences are available and can be included.
            - ✅ **Read counts included**  
              The OTU table contains the number of **sequence reads per OTU/ASV and sample**.
            - ✅ **Essential metadata available**  
              Metadata such as **sampling location, date, or environmental context** is provided.
            #### 🧹 Data Cleaning Requirements
            To ensure reliable downstream analyses, confirm the following:
            - ⚠️ **Controls are identified**  
              Non-sample entries (e.g., **negative controls, PCR blanks, positive controls**) are clearly identified and can be excluded, for example by removing them from the **Samples table**.
            - ⚠️ **Contaminants are removed**  
              Contaminant or unreliable detections (**OTUs/species/taxa**) are identified and can be excluded, for example by removing them from the **Taxonomy table**.
            """)

        if st.button("🔄 Convert Taxon Table to GBIF Upload", width='stretch'):
            if gbif_dataset_requirements == True:
                gbif_upload_conversion()
            else:
                st.info('Please read the GBIF dataset requirements first!')

        if st.button(label='🔗 GBIF Metabarcoding Data Toolkit (Test Environment)', width='stretch'):
            webbrowser.open('https://mdt.gbif-test.org/dataset/new')

        st.divider()

    ########################################################################################################################
    st.markdown("## 🇪🇺 Ecological Status Classes")
    with st.expander(expanded=False, label='See more'):
        st.success(
            "This section allows you to calculate widely used ecological status indices for freshwater ecosystems. "
            "You can compute the **European Fish Index (EFI)** and **Diathor Index**, as well as convert your taxon table into formats compatible with **Perlodes** and **Phylib** online tools. "
            "These indices help assess biodiversity, ecological quality, and ecosystem health based on species composition and abundance."
        )

        st.write('### European Fish Index')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            efi_measure = st.selectbox(label='Select Diversity Measure:', options=['Presence/Absence', 'Rel. Reads', f'{TTT_hash[:-1]} Diversity'], key='efi_measure')
        with col2:
            efi_placeholder = st.selectbox(label='Add Placeholder Data:', options=['Yes', 'No'], key='efi_placeholder')
        with col3:
            efi_info = st.toggle(label='Learn More About the EFI', key="efi_info")
        if efi_info == True:
            efi_table = path_to_ttt / 'WFD_conversion' / 'efi_table.xlsx'
            st.write('Table 1: Accepted variables for EFI calculation.')
            st.dataframe(pd.read_excel(efi_table, sheet_name='Sheet1').fillna(''), hide_index=True)
            st.write('Table 2: Available ecoregions for EFI calculation.')
            st.dataframe(pd.read_excel(efi_table, sheet_name='Sheet2').fillna(''), hide_index=True)
            efi_map = efi_table / 'WFD_conversion', 'efi_ecoregions.png'

        # --- Action Button ---
        if st.button(f"🔄 Calculate EFI", width='stretch'):
            calculate_efi_index()
        st.divider()

        st.write('### Diathor Index')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            diathor_measure = st.selectbox(label='Select Diversity Measure:', options=['Presence/Absence', 'Rel. Reads', f'{TTT_hash[:-1]} Diversity'], key='diathor_measure')
        with col2:
            pass
        with col3:
            pass
        # --- Action Button ---
        if st.button(f"🔄 Calculate Diathor Index", width='stretch'):
            calculate_diathor_index()
        st.divider()

        st.write('### Convert To Perlodes Format')
        # --- Initialize structure ---
        col1, col2 = st.columns([3,6])
        with col1:
            perlodes_measure = st.selectbox(label='Select Diversity Measure:', options=['Presence/Absence', 'Rel. Reads', f'{TTT_hash[:-1]} Diversity'], key='perlodes_measure')
        with col2:
            st.write('')
            st.write('')
            if st.button("🔗 Perlode Online Calculation", width='stretch'):
                pass
        # --- Action Button ---
        if st.button(f"🔄 Convert To Perlodes", width='stretch'):
            st.info('Coming soon!')
        st.divider()

        st.write('### Convert To Phylib Format')
        # --- Initialize structure ---
        col1, col2 = st.columns([3,6])
        with col1:
            phylib_measure = st.selectbox(label='Select Diversity Measure:', options=['Presence/Absence', 'Rel. Reads', f'{TTT_hash[:-1]} Diversity'], key='phylib_measure')
        with col2:
            st.write('')
            st.write('')
            if st.button("🔗 Phylib Online Calculation", width='stretch'):
                pass
        # --- Action Button ---
        if st.button(f"🔄 Convert To Phylib", width='stretch'):
            st.info('Coming soon!')
        st.divider()

else:
    st.write('')
    st.warning('Please Load a Taxon Table to continue!')