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
import sys
import subprocess
import webbrowser
import glob
from pathlib import Path
import plotly.express as px
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn2_unweighted
from matplotlib_venn import venn3
from matplotlib_venn import venn3_unweighted
from scipy.spatial import distance
import platform

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
    table_display_name = Path('/Volumes/Coruscant/test_folder/TTT_projects/eDNA_physalia_2026/TaXon_tables/eDNA_physalia_2026_taxon_table_merged_NCsub_fish.xlsx')
    df_taxon_table = pd.read_excel(table_display_name).fillna('')
    seq_loc = df_taxon_table.columns.tolist().index('Seq')
    hash_loc = df_taxon_table.columns.tolist().index('Hash')
    species_loc = df_taxon_table.columns.tolist().index('Species')
    df_taxon_cols = df_taxon_table.columns.tolist()[hash_loc+1:species_loc+1]
    df_samples = df_taxon_table.columns.tolist()[seq_loc+1:]
    trait_cols = ['Hash'] + df_taxon_table.columns.tolist()[species_loc+1:seq_loc]
    df_traits_table = df_taxon_table[trait_cols]
    taxon = df_taxon_table[df_taxon_cols]
    taxon_string = (df_taxon_table[df_taxon_cols].astype(str).agg(";".join, axis=1))
    similarity_loc = df_taxon_table.columns.get_loc("Similarity") + 1
    df_metadata_table = pd.read_excel(table_display_name, sheet_name='Metadata Table').fillna('')
    selected_metadata = 'Location2'

########################################################################################################################
# Page Setup
st.set_page_config(
    page_title="TaxonTableTools2",
    layout="wide",
    initial_sidebar_state="expanded"
)

########################################################################################################################
# general functions

def check_dependencies(tools=["cutadapt", "vsearch", "swarm", "blastn"]):
    missing = []
    for tool in tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        missing_tools = ', '.join(missing)
        if len(missing) == 1:
            st.error(f'WARNING: The following tool is missing: **{missing_tools}**')
        else:
            st.error(f'WARNING: The following tools are missing: **{missing_tools}**')
        st.warning('⚠️ Please install all required tools, either manually or using the "apscale_installer".')

def check_package_update(package_name):
    # Get the currently installed version
    try:
        installed_version = importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        print(f"{package_name} is not installed.")
        return

    # Check for updates
    res = update_check(package_name, installed_version)
    return res

def get_package_versions():
    packages = ["apscale", "apscale_gui", "apscale_blast", "boldigger3", "cutadapt"]

    for pkg in packages:
        try:
            version = importlib.metadata.version(pkg)
            st.write(f"{pkg}: {version}")
        except importlib.metadata.PackageNotFoundError:
            st.write(f"{pkg}: not installed")

    for tool in ['vsearch', 'swarm']:
        try:
            result = subprocess.run([tool, "--version"], capture_output=True, text=True, check=True)
            version = result.stderr.strip().split('\n')[0]
            st.write(f"{tool}: {version}")
        except:
            st.write(f"{tool}: not installed")

    try:
        result = subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=True)
        version = result.stdout.strip().split('\n')[0]
        st.write(version)
    except:
        st.write(f"blastn: not installed")

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
    table_name = table_display_name.stem
    active_project_path = st.session_state['active_project_path']
    os.makedirs(active_project_path / folder_name, exist_ok=True)

    if lib == 'matplot':
        file_pdf = active_project_path / folder_name / f'{table_name}_{suffix}.pdf'
        fig.savefig(file_pdf, dpi=300, bbox_inches="tight")
        st.success(f'Saved plot as pdf!')

    if lib == 'plotly':
        file_pdf = active_project_path / folder_name / f'{table_name}_{suffix}.pdf'
        file_html = active_project_path / folder_name / f'{table_name}_{suffix}.html'
        fig.update_layout(height=st.session_state['TTT_height'], width=st.session_state['TTT_width'])
        fig.write_image(file_pdf)
        fig.write_html(file_html)
        st.success(f'Saved plot as pdf and html!')

def export_table(folder_name, suffix, df, format):
    table_display_name = st.session_state['table_display_name']
    table_name = table_display_name.stem
    active_project_path = st.session_state['active_project_path']
    os.makedirs(active_project_path / folder_name, exist_ok=True)

    if format == 'xlsx':
        file_xlsx = active_project_path / folder_name / f'{table_name}_{suffix}.xlsx'
        df.to_excel(file_xlsx)
        st.success(f'Saved results as .xlsx!')

    if format == 'parquet':
        file_parquet = active_project_path / folder_name / f'{table_name}_{suffix}.parquet.snappy'
        df.to_parquet(file_parquet)
        st.success(f'Saved results as .parquet.snappy!')

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
        expected_columns = ["unique ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "E-value", "Flag", "Ambiguous taxa", "Status"]

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
        samples = merged_df.columns.tolist()[merged_df.columns.tolist().index('Seq') + 1:]
        metadata_df = pd.DataFrame([[i] + i.split('_') for i in samples])
        metadata_df.rename(columns={0: "Sample"}, inplace=True)

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
        metadata_df.rename(columns={0: "Sample"}, inplace=True)

        # Save table to a single Excel file
        to_excel_path = st.session_state['active_project_path'] / 'TaXon_tables' / str(
            st.session_state['taxon_table_creation_name'] + '_taxon_table.xlsx')
        with pd.ExcelWriter(to_excel_path, engine="openpyxl") as writer:
            merged_df.to_excel(writer, sheet_name="Taxon Table", index=False)
            metadata_df.to_excel(writer, sheet_name="Metadata Table", index=False)
        st.success(f'Saved Taxon Table as: {to_excel_path}')

########################################################################################################################
# Basic Stats

def basic_stats_reads():

    """ Calculate the number of reads per sample """
    samples = st.session_state['df_samples']
    df_taxon_table = st.session_state['df_taxon_table']
    y_values = {i:sum([j for j in df_taxon_table[i].values.tolist() if j != 0]) for i in samples}
    y_values = dict(sorted(y_values.items(), reverse=True, key=lambda item: item[1]))
    x_values = list(y_values.keys())

    max_reads, min_reads, avg_reads, stdev = max(y_values.values()), min(y_values.values()), round(statistics.mean(y_values.values()), 2), round(statistics.stdev(y_values.values()), 2)

    fig = go.Figure()
    fig = fig.add_trace(go.Bar(marker_color=st.session_state['TTT_color1'], x=x_values,y=list(y_values.values())))
    fig.update_yaxes(title='Reads', title_font=dict(size=st.session_state['TTT_fontsize']), tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_xaxes(title='Samples', title_font=dict(size=st.session_state['TTT_fontsize']), showticklabels=False, tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_layout(template=st.session_state['TTT_template'], font_size=st.session_state['TTT_fontsize'], title='Distribution of Reads')

    st.plotly_chart(fig, width='stretch', height='stretch')
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

    st.plotly_chart(fig, width='stretch')
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
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=[f'<i>{i}<i>' for i in list(species_occurrence_dict.keys())[:n]],  # Taxon names on x-axis
        y=list(species_occurrence_dict.values())[:n],
        marker_color=st.session_state['TTT_color1']
    ))

    fig.update_layout(template=st.session_state['TTT_template'], font_size=st.session_state['TTT_fontsize'], title=f'Top {n} Species')
    fig.update_yaxes(title = 'Occurrence', title_font=dict(size=st.session_state['TTT_fontsize']), tickfont=dict(size=st.session_state['TTT_fontsize']))
    fig.update_xaxes(title_font=dict(size=st.session_state['TTT_fontsize']), tickangle=-45, tickfont=dict(size=st.session_state['TTT_fontsize']*0.6))

    st.plotly_chart(fig, width='stretch', height='stretch')

    st.write(f'{len(all_species)} species across {len(samples)} samples')

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
    missing_replicates_handling = st.session_state['missing_replicates_handling']
    suffixes = st.session_state['suffixes']
    samples = st.session_state['df_samples']

    # --- derive clean sample names ---
    clean_samples = sorted({sample.rsplit('_', 1)[0] for sample in samples})

    columns_to_drop = []

    # --- merge replicates ---
    for sample in clean_samples:
        replicates = [f"{sample}_{suffix}" for suffix in suffixes]
        existing_replicates = [r for r in replicates if r in samples]
        # Case 1: all replicates present
        if len(existing_replicates) == len(suffixes):
            sub_df = df_taxon_table[existing_replicates]
            df_taxon_table[sample] = sub_df.apply(sum_if_less_than_X_zeros, axis=1, args=(cutoff_value,))
            columns_to_drop.extend(existing_replicates)
        # Case 2: missing replicates
        else:
            if missing_replicates_handling == "Keep":
                # do nothing → keep original replicate columns
                continue
            elif missing_replicates_handling == "Discard":
                columns_to_drop.extend(existing_replicates)

    # Drop replicate columns after loop
    df_taxon_table.drop(columns=columns_to_drop, inplace=True, errors="ignore")

    # Create new metadata df
    df_metadata_table = pd.DataFrame([[i, ''] for i in samples], columns=['Samples', 'Placeholder'])

    suffix = "merged"
    export_taxon_table(
        table_display_name,
        df_taxon_table,
        df_metadata_table,
        suffix
    )

    st.success(f'Replicates were merged and a new table was saved to "{str(table_display_name.name).replace(".xlsx", "_merged.xlsx")}"')

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

        export_plot('Venn_diagrams', f'{venn_metadata}_{venn_taxon_level}', plt, 'matplot')

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
    all_taxa = df_taxon_table_simple[fig_taxon].unique()

    # Output Table
    # TBD

    # Output Plot
    fig = go.Figure()
    colors = get_colors_from_scale(st.session_state['TTT_colorscale'], len(all_taxa))
    for c, taxon in enumerate(all_taxa):
        y_values = df_taxon_table_simple_s[df_taxon_table_simple_s[fig_taxon] == taxon][df_samples].values.tolist()[0]
        x_values = df_samples
        fig.add_trace(go.Bar(x=x_values, y=y_values, name=taxon, marker=dict(color=colors[c])))

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
    st.plotly_chart(fig, width='stretch', height='stretch')

    # Export
    export_plot("Read_proportions_plots", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", fig, "plotly")

########################################################################################################################
# Alpha Diversity

def alpha_richness_diagram():

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
    fig_color = st.session_state['richness_color'] # color sequence or scale
    fig_mode = st.session_state['richness_mode'] # Bar or Box/Violin
    fig_level = st.session_state['richness_level'] # sample or metadata
    fig_taxon = st.session_state['richness_taxon'] # taxonomic level
    fig_nan = st.session_state['richness_nan'] # include or exclude
    fig_boxpoints = st.session_state['richness_boxpoints']

    if fig_mode == 'Bar':

        # Handle Metadata
        if fig_level != 'Samples':
            df_taxon_table, df_samples = concatenate_by_metadata(fig_level)
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)

        # Analysis
        if fig_nan == 'Exclude':
            df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']
        df_taxon_table_simple_s = (df_taxon_table_simple[[fig_taxon] + df_samples].groupby(fig_taxon, as_index=False).sum())

        # Output Table
        # TBD

        fig = go.Figure()
        if fig_color == 'Single Color':
            colors = [st.session_state['TTT_color1'] ] * len(df_samples)
        else:
            colors = get_colors_from_scale(st.session_state['TTT_colorscale'], len(df_samples))
        for c, sample in enumerate(df_samples):
            x_values = [sample]
            y_values = [len([i for i,j in df_taxon_table_simple_s[[fig_taxon, sample]].values.tolist() if j != 0])]
            fig.add_trace(go.Bar(x=x_values, y=y_values, name=sample, marker=dict(color=colors[c])))

        # Specific Layout
        if fig_taxon == 'Hash':
            fig_taxon = st.session_state['TTT_hash'][:-1]
        fig.update_xaxes(dtick='linear')
        fig.update_yaxes(title=f'{fig_taxon} Richness')
        fig.update_layout(title=f'{fig_taxon} Richness')

        # Default Layout
        fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)

        # Plot
        st.plotly_chart(fig)

        # Export
        export_plot("Alpha_diversity", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", fig, "plotly")

    elif fig_mode == 'Box' or fig_mode == 'Violin':
        if fig_level == 'Samples':
            st.error('Please select a metadata category instead of "Samples"!')
            return

        # Collect metadata
        test_metadata = [i for i in list(df_metadata_table[fig_level].unique()) if i != '']

        # Create simple table
        df_taxon_table_simple = simple_taxon_table(df_taxon_table)

        # Analysis
        if fig_nan == 'Exclude':
            df_taxon_table_simple = df_taxon_table_simple[df_taxon_table_simple[fig_taxon] != '']
        df_taxon_table_simple_s = (df_taxon_table_simple[[fig_taxon] + df_samples].groupby(fig_taxon, as_index=False).sum())

        # Output Table
        # TBD

        fig = go.Figure()

        if fig_color == 'Single Color':
            colors = [st.session_state['TTT_color1']] * len(test_metadata)
        else:
            colors = get_colors_from_scale(st.session_state['TTT_colorscale'], len(test_metadata))

        for c, metadata in enumerate(test_metadata):
            test_samples = list(df_metadata_table[df_metadata_table[fig_level] == metadata]['Samples'].unique())
            if len(test_samples) == 1:
                test_samples = [test_samples]
            y_values = [(df_taxon_table_simple_s[sample] != 0).sum() for sample in test_samples]
            if fig_mode == 'Box':
                fig.add_trace(go.Box(y=y_values, name=metadata, marker=dict(color=colors[c])))
            if fig_mode == 'Violin':
                fig.add_trace(go.Violin(y=y_values, name=metadata, marker=dict(color=colors[c])))

        # Specific Layout
        if fig_taxon == 'Hash':
            fig_taxon = st.session_state['TTT_hash'][:-1]
        fig.update_xaxes(dtick='linear')
        fig.update_yaxes(title=f'{fig_taxon} Richness', rangemode='tozero')
        fig.update_layout(title=f'{fig_taxon} Richness')
        if fig_boxpoints == 'All' and fig_mode == 'Box':
            fig.update_traces(boxpoints='all', jitter=st.session_state['TTT_jitter'])
        if fig_boxpoints == 'All' and fig_mode == 'Violin':
            fig.update_traces(points='all', jitter=st.session_state['TTT_jitter'])

        # Default Layout
        fig.update_yaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_xaxes(title_font=dict(size=fontsize), tickfont=dict(size=fontsize))
        fig.update_layout(template=template, font_size=fontsize, showlegend=showlegend)

        # Plot
        st.plotly_chart(fig)

        # Export
        export_plot("Alpha_diversity", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", fig, "plotly")
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
    df_taxon_table_simple = simple_taxon_table(df_taxon_table)

    # Options
    df_metadata_table = st.session_state['df_metadata_table'].copy()
    df_samples = st.session_state['df_samples'].copy()
    fig_labels = st.session_state['betadiv_labels'] # color sequence or scale
    fig_mode = st.session_state['betadiv_mode'] # spearman or bray curtis
    fig_measure = st.session_state['betadiv_measure'] # spearman or bray curtis
    fig_level = st.session_state['betadiv_level'] # sample or metadata
    fig_taxon = st.session_state['betadiv_taxon'] # taxonomic level
    fig_nan = st.session_state['betadiv_nan'] # include or exclude

    # Handle Metadata
    if fig_level != 'Samples':
        df_taxon_table, df_samples = concatenate_by_metadata(fig_level)
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
    st.plotly_chart(fig)

    # Export
    export_plot("Beta_diversity", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", fig, "plotly")
    export_table("Beta_diversity", f"{name}_{fig_mode}_{fig_level}_{fig_taxon}", df, "xlsx")

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

    available_templates_list = ['seaborn', 'ggplot2', 'simple_white', 'plotly', 'plotly_dark', 'presentation', 'plotly_white']

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
    TTT_jitter = user_data_df[user_data_df['Variable'] == 'TTT_jitter']['Value'].values.tolist()[0]
    TTT_color1 = user_data_df[user_data_df['Variable'] == 'TTT_color1']['Value'].values.tolist()[0]
    TTT_color2 = user_data_df[user_data_df['Variable'] == 'TTT_color2']['Value'].values.tolist()[0]
    TTT_color3 = user_data_df[user_data_df['Variable'] == 'TTT_color3']['Value'].values.tolist()[0]
    TTT_colorsequence = user_data_df[user_data_df['Variable'] == 'TTT_colorsequence']['Value'].values.tolist()[0]
    TTT_colorscale = user_data_df[user_data_df['Variable'] == 'TTT_colorscale']['Value'].values.tolist()[0]

########################################################################################################################

with st.sidebar:
    #####################################
    st.write(""" # Projects """)

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
        if st.button('Open Active Project', width='stretch'):
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
            if st.button('Load Table', width='stretch'):
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
                st.session_state['df_samples'] = st.session_state['df_taxon_table'].columns.tolist()[seq_loc+1:]
                species_loc = st.session_state['df_taxon_table'].columns.tolist().index('Species')
                # Traits -> df_traits_table
                st.session_state['df_traits_table'] = st.session_state['df_taxon_table'].columns.tolist()[species_loc+1:seq_loc]
                # Hashes -> df_hashes
                st.session_state['df_hashes'] = st.session_state['df_taxon_table']['Hash']
                # Taxon columns
                st.session_state['df_taxon_cols'] = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

            if 'df_taxon_table' in st.session_state and 'table_display_name' in st.session_state:
                st.success(f'Active table: {st.session_state["table_display_name"].stem}')
                if st.button('Open Active Table', width='stretch'):
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
            # SCATTER SIZE
            st.session_state['TTT_jitter'] = st.sidebar.number_input('Jitter', 0.0, 1.0, 0.3)

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

            if st.button("🔄 Save Settings", width='stretch'):
                settings_df.to_excel(user_preferences_xlsx, index=False)
                st.success('Settings were saved!')


########################################################################################################################
st.markdown("## 🧬 Create Taxon Table")
if 'df_taxon_table' not in st.session_state or 'table_display_name' not in st.session_state:
    expanded=True
else:
    expanded=False

if path_to_projects.exists():

    with st.expander(expanded=expanded, label='See more'):
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
            basic_stats_reads()
            basic_stats_OTUs()
            top_n_taxa_plot()

    ########################################################################################################################
    st.markdown("## 🛠️ Table Processing")
    with st.expander(expanded=False, label='See more'):

        st.write('### Replicate Merging')
        # --- Derive sample structure ---
        df_samples = st.session_state['df_samples']
        hash_label = st.session_state['TTT_hash']
        suffixes = sorted({i.split('_')[-1] for i in df_samples})
        st.session_state['suffixes'] = suffixes
        unique_samples = sorted({'_'.join(i.split('_')[:-1]) for i in df_samples})
        st.session_state['unique_samples'] = unique_samples
        expected_samples = [f"{s}_{suf}" for suf in suffixes for s in unique_samples]
        n_suffixes = len(suffixes)
        cutoff_options = list(range(1, n_suffixes + 1))[::-1]
        # --- Overview Metrics ---
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("Total Samples", len(df_samples))
        col2.metric("Unique Samples", len(unique_samples))
        col3.metric("Replicate Suffixes", ", ".join(suffixes))
        # --- Cutoff Selection ---
        with col4:
            cutoff_value = st.selectbox(label=f"{hash_label} must be present in at least:",options=cutoff_options, format_func=lambda x: f"{x} replicate(s)", key='cutoff_value')
        with col5:
            missing_replicates_handling = st.selectbox(label=f"Handle samples which miss replicates:",options=['Keep', 'Remove'], key='missing_replicates_handling')
        if cutoff_value == 1:
            st.info(f"All {hash_label} will be kept after replicate merging.")
        else:
            st.warning(f"Only {hash_label} present in "f"{cutoff_value}/{n_suffixes} replicates "f"will be kept after replicate merging.")
        # --- Action Button ---
        if st.button("🔄 Merge Replicates", use_container_width=True):
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
        if st.button("🔄 Subtract Negative Controls", use_container_width=True):
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
        if st.button("🔄 Apply read-based filter", use_container_width=True):
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
        if st.button("🔄 Apply read-based normalisation", use_container_width=True):
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
        if st.button("🔄 Apply taxonomy filter", use_container_width=True):
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
        if st.button("🔄 Apply sample filter", use_container_width=True):
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
        if st.button("🔄 Apply trait filter", use_container_width=True):
            trait_filter()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🔁 Table Conversion")
    with st.expander(expanded=False, label='See more'):

        st.write('### Simplify Taxon Table')
        # --- Initialize structure ---
        st.info(
            f"{st.session_state['TTT_hash']} with identical taxonomic assignments "
            f"will be merged into a single entry. The representative "
            f"{st.session_state['TTT_hash']} will be selected based on the highest "
            f"read abundance and sequence similarity."
        )
        # --- Action Button ---
        if st.button("🔄 Simplify Taxon Table", use_container_width=True):
            simplify_table()
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
        if st.button("🔄 Add Traits From File", use_container_width=True):
            add_traits_from_file()
        st.divider()


        st.write('### Sort Samples')
        # --- Initialize structure ---
        st.info('Sort your Metadata Table as desired. The samples of the Taxon Table will then be sorted accordingly!')
        # --- Action Button ---
        if st.button("🔄 Sort Samples", use_container_width=True):
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
        if st.button("🔄 Rename Samples", use_container_width=True):
            sample_rename()
        st.divider()


        st.write('### Export Table To Fasta')
        # --- Initialize structure ---
        st.info(f"{st.session_state['TTT_hash']} will be exported to .FASTA format.")
        # --- Action Button ---
        if st.button("🔄 Export Table To Fasta", use_container_width=True):
            export_fasta()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🧪 Sample Comparison")
    with st.expander(expanded=False, label='See more'):

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
        elif st.button(f"🔄 Generate Venn{len(venn_categories)} Diagram", use_container_width=True):
            venn_diagram()
        st.divider()

    ########################################################################################################################
    st.markdown("## 📈 Read Proportions")
    with st.expander(expanded=False, label='See more'):

        st.write('### Read Distribution Diagram')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            readdist_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols'], index=2, key='readdist_taxon')
            readdist_level = st.selectbox(label='Select Sample Level:', options=['Samples'] + available_metadata, key='readdist_level')
        with col2:
            readdist_mode = st.selectbox(label='Select Display Mode:', options=['Relative Reads', 'Absolute Reads'], key='readdist_mode')
            readdist_nan = st.selectbox(label='Handle Missing Taxonomy:', options=['Exclude', 'Include'], key='readdist_nan')
        with col3:
            readdist_color = st.selectbox(label='Select Color Mode:', options=['Color scale'], key='readdist_color')
        # --- Action Button ---
        if st.button(f"🔄 Generate Read Distribution Diagram", use_container_width=True):
            readdist_diagram()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🌱 Alpha Diversity")
    with st.expander(expanded=False, label='See more'):

        st.write('### Taxon Richness Diagram')
        # --- Initialize structure ---
        col1, col2, col3 = st.columns(3)
        with col1:
            richness_taxon = st.selectbox(label='Select Taxonomic Level:', options=st.session_state['taxa_cols2'], index=6, key='richness_taxon')
            richness_level = st.selectbox(label='Select Sample Level:', options=['Samples'] + available_metadata, key='richness_level')
        with col2:
            richness_mode = st.selectbox(label='Select Display Mode:', options=['Bar', 'Box', 'Violin'], key='richness_mode')
            richness_boxpoints = st.selectbox(label='Boxpoints Display:', options=['All', 'None'], key='richness_boxpoints')
        with col3:
            richness_nan = st.selectbox(label='Handle Missing Taxonomy:', options=['Exclude', 'Include'], key='richness_nan')
            richness_color = st.selectbox(label='Select Color Mode:', options=['Color scale', 'Single Color'], key='richness_color')
        # --- Action Button ---
        if st.button(f"🔄 Generate {richness_taxon} Richness Diagram", use_container_width=True):
            alpha_richness_diagram()
        st.divider()

    ########################################################################################################################
    st.markdown("## 🌎 Beta Diversity")
    with st.expander(expanded=False, label='See more'):

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
        if st.button(f"🔄 Generate {betadiv_calculation} {betadiv_measure} Diagram", use_container_width=True):
            beta_diversity_matrix()
        st.divider()

    ########################################################################################################################
    st.markdown("## 📆 Time Series")
    with st.expander(expanded=False, label='See more'):
        st.info('Coming soon!')

    ########################################################################################################################
    st.markdown("## 🔗 API Modules")
    with st.expander(expanded=False, label='See more'):
        st.info('Coming soon!')

    ########################################################################################################################
    st.markdown("## 🇪🇺 Ecological Status Classes")
    with st.expander(expanded=False, label='See more'):
        st.info('Coming soon!')
else:
    st.write('')
    st.warning('Please Load a Taxon Table to continue!')