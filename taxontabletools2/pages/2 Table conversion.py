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
from taxontabletools2.table_conversion import presence_absence_table
from taxontabletools2.table_conversion import simplify_table
from taxontabletools2.table_conversion import add_traits_from_file
from taxontabletools2.table_conversion import add_traits_from_file
from taxontabletools2.table_conversion import rename_samples
from taxontabletools2.table_conversion import sort_samples
from taxontabletools2.table_conversion import merge_ESV_tables
from taxontabletools2.table_conversion import open_reading_frame
import streamlit as st
from taxontabletools2.start import start
from Bio.Seq import CodonTable

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## Table conversion """)

# Ensure the session state variables are retained across subpages
if 'TaXon_table_df' not in st.session_state:
    error1, error2, error3 = st.columns(3)
    with error2:
        st.write('##### Please load a TaXon table to continue!')

elif 'TaXon_table_xlsx' not in st.session_state:
    error1, error2, error3 = st.columns(3)
    with error2:
        st.write('##### Please load a TaXon table to continue!')

elif st.session_state['TaXon_table_df'].empty == True:
    error1, error2, error3 = st.columns(3)
    with error2:
        st.write('##### Please check your TaXon table. It seems to be empty!')

elif st.session_state['TaXon_table_xlsx'].is_file() == False:
    error1, error2, error3 = st.columns(3)
    with error2:
        st.write('##### The TaXon table does not exist in the selected project folder!')

########################################################################################################################

else:
    ########################################################################################################################
    st.markdown("---")
    st.write('### Presence absence table')
    with st.container():
        b1, b2 = st.columns(2)
        with b1:
            st.write('Converts reads to presence/absence format (1/0).')
        with b2:
            pass

        if st.button('Convert to presence/absence table'):
            presence_absence_table(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'])
            st.success('Converted table to presence/absence data. Please remember to reload the TaXon table!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### Simplify table')
    with st.container():
        d1, d2 = st.columns(2)
        with d1:
            st.write(f'This module will merge all {st.session_state["clustering_unit"]} with identical taxonomic assignments into a single operational taxonomic unit, called \"spOTU\".')
            st.write('Please not that previously assigned traits must be re-added again. Also, sequences will be removed and the similarity will be set to 100.')
        with d2:
            pass
        if st.button('Simplify table'):
            simplify_table(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'])
            st.success('Simplified table. Please remember to reload the TaXon table!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### Add traits from file')
    with st.container():
        e1, e2 = st.columns(2)
        with e1:
            new_traits_file = st.file_uploader('Upload your .xlsx file that contains the traits!')
        with e2:
            if new_traits_file != None:
                st.session_state['new_traits_df'] = pd.read_excel(new_traits_file)
                st.selectbox('Select taxonomic level:',
                             [i for i in st.session_state['new_traits_df'] if i in available_taxonomic_levels_list])
                st.multiselect('Select traits to add:',
                             [i for i in st.session_state['new_traits_df'] if i not in available_taxonomic_levels_list])

        if st.button('Add traits to table'):
            pass
            st.success('Text')

    ########################################################################################################################
    st.markdown("---")
    st.write('### Sort samples')
    with st.container():
        st.write('Sort samples based on the current sorting of the metadata table.')
        if st.button('Sort samples'):
            sort_samples(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'])
            st.success('Samples were sorted according to the metadata table!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### Rename samples')
    with st.container():
        g1, g2 = st.columns(2)
        with g1:
            available_metadata = st.session_state['metadata_df'].columns.tolist()[1:]
            selected_metadata = st.selectbox('Select metadata to test:', available_metadata)
        with g2:
            pass

        st.write(f'Rename samples based on the "{selected_metadata}" column.')
        if st.button('Rename samples'):
            rename_samples(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'], selected_metadata)
            st.success(f'Samples were renamed based on the "{selected_metadata}" column!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### Merge ESV tables')
    with st.container():
        h1, h2 = st.columns(2)
        with h1:
            ESV_table_2 = st.file_uploader('Load another TaXon table (must contain unique ESVs)')
        with h2:
            suffix = st.text_input('Enter a suffix for the merged table:')

        if st.button('Merge ESV tables'):
            # Load df2
            df2 = pd.read_excel(ESV_table_2).fillna('')
            df2_metadata = pd.read_excel(ESV_table_2, sheet_name='Metadata Table').fillna('')

            merge_ESV_tables(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], df2, st.session_state['metadata_df'], df2_metadata, suffix)

    ########################################################################################################################
    st.markdown("---")
    st.write('### Aminoc acid translation')
    with st.container():
        k1, k2 = st.columns(2)
        with k1:
            all_tables = {i[1].names[0]:i[1].id for i in CodonTable.unambiguous_dna_by_id.items()}
            selection = st.selectbox('Select translation table:', list(all_tables.keys()))
        with k2:
            orf_filter = st.selectbox('Filter sequences for open reading frame:', ['No', 'Yes'])

        if st.button('Translate sequences'):
            table_number = all_tables[selection]
            open_reading_frame(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['metadata_df'], st.session_state['traits_df'], table_number, orf_filter)


