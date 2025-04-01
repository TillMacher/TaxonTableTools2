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
import streamlit as st
from taxontabletools2.start import start

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## Table processing """)

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
    st.write('### 1) Replicate merging ')
    with st.container():
        replicate_dict = collect_replicates(st.session_state['samples'])

        n_samples = len(replicate_dict[0])
        n_replicates_suffixes = len(replicate_dict[1])
        a1, a2 = st.columns(2)
        with a1:
            replicate_selector = st.selectbox('Replicates select:', ['Manual', 'Detect from file'])

            if replicate_selector == 'Detect from file':
                st.write(f'Replicate suffixes: {replicate_dict[1]}')
            else:
                replicates_suffixes = st.text_input('Enter replicates suffixes (separate by comma):', 'A,B')
                replicates_suffixes = [i for i in replicates_suffixes.split(',')]
                st.write(f'Replicate suffixes: {replicates_suffixes}')
                n_replicates_suffixes = len(replicates_suffixes)
                replicate_dict[1] = replicates_suffixes

        with a2:
                replicate_cutoff = st.selectbox('Replicates cutoff:', [i for i in range(1, n_replicates_suffixes + 1)][::-1])
                st.write(f'OTU must be present in {replicate_cutoff} of {n_replicates_suffixes} replicates!')

                replicate_non_matching_samples = st.selectbox('Samples without matching replicates:', ['Keep (without merging)', 'Remove from dataset'])

        test, missing_samples = check_replicates(replicate_dict, st.session_state['samples'])
        if test == False:
            st.error('Warning: Not all samples are present with all expected replicates. Please check you table!')
            st.error(', '.join(list(missing_samples)[:20]))

        if st.button('Merge replicates'):
            replicate_merging(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'], replicate_cutoff, replicate_dict, replicate_non_matching_samples)
            st.success('Merged replicates. Remember to reload the TaXon table!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### 2) Negative control subtraction ')
    with st.container():
        b1, b2 = st.columns(2)
        with b1:
            search_mask = st.text_input('Enter negative control identifier:', 'NC_')
        with b2:
            mode = st.selectbox('Select how to deal with negatice control subtraction', ['Sum of all NCs', 'Maximum in all NCs', 'Average in all NCs',])
        samples_to_remove = st.multiselect('Selected NCs:', [i for i in st.session_state['samples'] if search_mask in i], default=[i for i in st.session_state['samples'] if search_mask in i])
        if st.button('Subtract controls'):
            negative_control_subtraction(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'], mode, samples_to_remove)
            st.success(f'Subtracted reads from {len(samples_to_remove)} negative controls. Remember to reload the TaXon table!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### 3) Read-based filter ')
    with st.container():
        c1, c2 = st.columns(2)
        with c1:
            read_filter_mode = st.selectbox('Read-based filter mode:', ['Relative', 'Absolute'])
            if read_filter_mode == 'Absolute':
                read_filter_threshold = st.number_input('Cutoff:', min_value=0, value=1000)
            else:
                read_filter_threshold = st.number_input('Cutoff (%):', min_value=0.000, max_value=100.0, value=0.01)
            if st.button('Read-based filter'):
                read_based_filter(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'], read_filter_mode, read_filter_threshold)
                st.success(f'Filtered samples according to reads. Remember to reload the TaXon table!')
        with c2:
            read_normalization = st.selectbox('Normalization filter mode:', ['Average', 'Minimum', 'Custom'])
            if read_normalization == 'Average':
                n_reads = round(st.session_state['TaXon_table_df'][st.session_state['samples']].sum().mean())
                normalization_threshold = st.number_input('Normalization read threshold:', min_value=0, value=n_reads)
            elif read_normalization == 'Minimum':
                n_reads = min(st.session_state['TaXon_table_df'][st.session_state['samples']].sum())
                normalization_threshold = st.number_input('Normalization read threshold:', min_value=0, value=n_reads)
            else:
                normalization_threshold = st.number_input('Normalization read threshold:', min_value=0,
                                                          value=60000)
            if st.button('Read-based normalization'):
                read_based_normalisation(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'], normalization_threshold)
                st.success(f'Normalized samples according to reads. Remember to reload the TaXon table!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### 4) Taxonomic filtering ')
    with st.container():
        d1, d2 = st.columns(2)
        with d1:
            taxonomic_level = st.selectbox('Select a taxonomic level:', st.session_state['available_taxonomic_levels_list'])
            available_taxa = sorted(set(st.session_state['TaXon_table_df'][taxonomic_level].values.tolist()))
        with d2:
            selector = st.selectbox('Select mode:', ['Select all available taxa', 'Select no taxa'])
        if selector == 'Select all available taxa':
            taxa_to_remove = st.multiselect(' ', available_taxa, default=available_taxa)
        else:
            taxa_to_remove = st.multiselect(' ', available_taxa)
        taxonomic_filtering_suffix = st.text_input('Enter suffix for the filtered table:', '', key='taxonomic_filtering_suffix')
        if st.button('Remove selected taxa'):
            mode= 'Remove'
            taxonomic_filtering(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'],
                                st.session_state['samples'], st.session_state['metadata_df'],
                                st.session_state['traits_df'], mode, taxa_to_remove, taxonomic_level, taxonomic_filtering_suffix)
            st.success(f'Removed {len(taxa_to_remove)} taxa. Remember to reload the TaXon table!')
        if st.button('Keep selected taxa'):
            mode = 'Keep'
            taxonomic_filtering(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'],
                                st.session_state['samples'], st.session_state['metadata_df'],
                                st.session_state['traits_df'], mode, taxa_to_remove, taxonomic_level, taxonomic_filtering_suffix)
            st.success(f'Kept {len(taxa_to_remove)} taxa. Remember to reload the TaXon table!')

    ########################################################################################################################
    st.markdown("---")
    st.write('### 5) Sample filtering ')
    with st.container():
        e1, e2 = st.columns(2)
        with e1:
            search_mask = st.text_input('Type to filter for samples:', '', key='sample_filtering_search_mask')
        with e2:
            selector = st.selectbox('Select mode:', ['Select all available samples', 'Select no samples'])
        samples_to_remove = st.multiselect('Selected samples:', [i for i in st.session_state['samples'] if search_mask in i], default=[i for i in st.session_state['samples'] if search_mask in i])
        sample_filtering_suffix = st.text_input('Enter suffix for the filtered table:', '', key='sample_filtering_suffix')
        if st.button('Remove selected samples'):
            mode = 'Remove'
            sample_filtering(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'],
                                st.session_state['samples'], st.session_state['metadata_df'],
                                st.session_state['traits_df'], mode, samples_to_remove, sample_filtering_suffix)
            st.success(f'Removed {len(samples_to_remove)} samples. Remember to reload the TaXon table!')
        if st.button('Keep selected samples'):
            mode = 'Keep'
            sample_filtering(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'],
                                st.session_state['samples'], st.session_state['metadata_df'],
                                st.session_state['traits_df'], mode, samples_to_remove, sample_filtering_suffix)
            st.success(f'Kept {len(samples_to_remove)} samples. Remember to reload the TaXon table!')