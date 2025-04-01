import streamlit as st
import pandas as pd
import numpy as np
import subprocess
import glob, sys, time, statistics, os.path, random
from pathlib import Path
import psutil
import plotly.graph_objects as go
from taxontabletools2.start import TTT_variables
from taxontabletools2.utilities import load_df
from taxontabletools2.utilities import collect_replicates
from taxontabletools2.utilities import collect_sample_stats
from taxontabletools2.basic_stats import basic_stats_reads
from taxontabletools2.basic_stats import basic_stats_OTUs
from taxontabletools2.basic_stats import taxonomic_resolution
from taxontabletools2.basic_stats import taxonomic_richness
from taxontabletools2.utilities import collect_traits
from taxontabletools2.utilities import strip_traits
from taxontabletools2.utilities import check_replicates
from taxontabletools2.table_processing import replicate_merging
from taxontabletools2.table_processing import negative_control_subtraction
from taxontabletools2.table_processing import read_based_filter
from taxontabletools2.table_processing import read_based_normalisation
from taxontabletools2.table_processing import taxonomic_filtering
from taxontabletools2.table_processing import sample_filtering
from taxontabletools2.start import start

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## Taxon Table conversion """)

def merge_otu_and_taxonomy_tables():
    """
    Merges an OTU table and a taxonomy table based on unique IDs.
    Ensures proper formatting, column alignment, and metadata creation.
    """

    # Load OTU table
    otu_table = pd.read_excel(st.session_state['otu_table_path']).fillna('')

    # Determine the correct sheet name for the taxonomy table
    taxonomy_sheets = {'APSCALE': 'Taxonomy table', 'BOLDigger': 'Sheet1'}
    sheet_name = taxonomy_sheets.get(st.session_state['taxonomy_table_format'], None)

    if sheet_name is None:
        st.error('Invalid taxonomy table format.')
        return

    # Load taxonomy table
    try:
        taxonomy_table = pd.read_excel(st.session_state['taxonomy_table_path'], sheet_name=sheet_name).fillna('')
    except ValueError:
        st.error(f"Worksheet named '{sheet_name}' not found")
        return

    # Extract unique IDs from both tables
    otu_table_IDs = otu_table['unique_ID'].tolist()
    taxonomy_column_map = {'APSCALE': 'unique ID', 'BOLDigger': 'id'}
    taxonomy_id_col = taxonomy_column_map.get(st.session_state['taxonomy_table_format'])

    if taxonomy_id_col in taxonomy_table.columns:
        taxonomy_table_IDs = taxonomy_table[taxonomy_id_col].tolist()
        taxonomy_table.rename(columns={taxonomy_id_col: 'unique_ID'}, inplace=True)

        # Additional renaming for BOLDigger format
        if st.session_state['taxonomy_table_format'] == 'BOLDigger':
            taxonomy_table.rename(columns={'pct_identity': 'Similarity'}, inplace=True)

        cont = True
    else:
        cont = False

    # Check for ID mismatches
    if otu_table_IDs != taxonomy_table_IDs:
        st.error(
            'WARNING: The IDs of the Read table and Taxonomy table do not match! Make sure they are identically sorted!')
        return
    elif not cont:
        st.error('WARNING: Unable to find the required sheets in the provided taxonomy table.')
        return

    # Merge OTU and taxonomy tables on unique_ID
    taxon_table = taxonomy_table.merge(otu_table, how='inner', on='unique_ID')

    # Move 'Seq' column to correct position
    first_sample = otu_table.columns[2]
    loc_first_sample = taxon_table.columns.get_loc(first_sample)
    seq = taxon_table.pop('Seq')  # Remove 'Seq' column
    taxon_table.insert(loc_first_sample, 'Seq', seq)  # Insert at new position

    # Rename 'unique_ID' to 'ID'
    taxon_table.rename(columns={'unique_ID': 'ID'}, inplace=True)

    # Create metadata DataFrame
    samples = otu_table.columns[2:-1]
    metadata_df = pd.DataFrame({'Sample': samples, 'Metadata': ''})

    # Ensure 'Kingdom' column exists
    if 'Kingdom' not in taxon_table.columns:
        loc_phylum = taxon_table.columns.get_loc('Phylum')
        taxon_table.insert(loc_phylum, 'Kingdom', 'Life')

    # Merge 'Genus' and 'Species' columns if required
    if st.session_state['genus_species_merge'] == 'Merge':
        taxon_table['Species'] = taxon_table.apply(
            lambda row: f"{row['Genus']} {row['Species']}" if row['Genus'] and row['Species'] else row['Species'],
            axis=1
        )

    # Display merged table
    st.success('Successfully merged your tables:')
    st.dataframe(taxon_table)

    # Save merged table and metadata to an Excel file
    with pd.ExcelWriter(st.session_state['taxon_table_path']) as writer:
        taxon_table.to_excel(writer, sheet_name='Taxon Table', index=False)
        metadata_df.to_excel(writer, sheet_name='Metadata Table', index=False)

    st.success(f'Saved as "{st.session_state["taxon_table_path"]}"')

a1, a2 = st.columns(2)
with a1:
    st.session_state['otu_table_path'] = st.file_uploader('Select your Read table:')
    st.session_state['otu_table_format'] = st.selectbox('Read table format', ['APSCALE'])

    if st.session_state['otu_table_path'] == None:
        st.write('Please select a Read table!')

with a2:
    st.session_state['taxonomy_table_path'] = st.file_uploader('Select your Taxonomy table:')
    st.session_state['taxonomy_table_format'] = st.selectbox('Taxonomy table format', ['APSCALE', 'BOLDigger'])
    st.session_state['genus_species_merge'] = st.selectbox('Merge the genus and species epithet column', ['Merge', 'No merge'], index=1)

    if st.session_state['taxonomy_table_path'] == None:
        st.write('Please select a Taxonomy table!')

## Save table to project folder
st.session_state['taxon_table_name'] = st.text_input('Taxon table name',
                                                     f'{st.session_state["path_to_outdirs"].name}_taxon_table') + '.xlsx'
st.session_state['taxon_table_path'] = st.session_state['path_to_outdirs'].joinpath('TaXon_tables',
                                                                                  st.session_state['taxon_table_name'])

if st.session_state['taxonomy_table_path'] != None and st.session_state['otu_table_path'] != None:
    if st.button('Convert to TaXon table format'):
        merge_otu_and_taxonomy_tables()















