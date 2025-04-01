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
from taxontabletools2.utilities import export_plot
from taxontabletools2.table_processing import replicate_merging
from taxontabletools2.table_processing import negative_control_subtraction
from taxontabletools2.table_processing import read_based_filter
from taxontabletools2.table_processing import read_based_normalisation
from taxontabletools2.table_processing import taxonomic_filtering
from taxontabletools2.table_processing import sample_filtering
from taxontabletools2.sample_comparison import venn
from taxontabletools2.read_analyses import read_props_bar
from taxontabletools2.read_analyses import read_props_heatmap
import streamlit as st
from taxontabletools2.start import start

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## Read proportions """)

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
    st.write('### Read proportion analyses')

    st.write('### Bar chart')

    reads1, reads2 = st.columns(2)

    with reads1:
        available_metadata = ['All samples']  + st.session_state['metadata_df'].columns.tolist()[1:]
        selected_metadata = st.selectbox('Select metadata to test:', available_metadata)

    with reads2:
        taxonomic_level = st.selectbox('Select taxonomic level:', st.session_state['available_taxonomic_levels_list'], index=5)

    if st.button(f'Calculate {taxonomic_level} reads proportions for {selected_metadata}', key='readprops_bar'):

        fig_read_props_1 = read_props_bar(
            st.session_state['path_to_outdirs'],
            st.session_state['TaXon_table_xlsx'],
            st.session_state['TaXon_table_df'],
            st.session_state['samples'],
            st.session_state['metadata_df'],
            selected_metadata,
            st.session_state['traits_df'],
            '',
            users_settings,
            {'selected_metadata': selected_metadata, 'taxonomic_level': taxonomic_level}
        )

        st.plotly_chart(fig_read_props_1, use_container_width=True, config=st.session_state['config'])
        export_plot(st.session_state['TaXon_table_xlsx'], 'Read_proportions_plots', 'bar_chart', selected_metadata, fig_read_props_1, users_settings, 'plotly')

    st.write('### Heatmap')

    reads3, reads4 = st.columns(2)

    with reads3:
        available_metadata = ['All samples']  + st.session_state['metadata_df'].columns.tolist()[1:]
        selected_metadata = st.selectbox('Select metadata to test: ', available_metadata)

    with reads4:
        taxonomic_level = st.selectbox('Select taxonomic level: ', st.session_state['available_taxonomic_levels_list'], index=5)

    if st.button(f'Calculate {taxonomic_level} reads proportions for {selected_metadata}', key='readprops_heatmap'):

        fig_read_props_1 = read_props_heatmap(
            st.session_state['path_to_outdirs'],
            st.session_state['TaXon_table_xlsx'],
            st.session_state['TaXon_table_df'],
            st.session_state['samples'],
            st.session_state['metadata_df'],
            selected_metadata,
            st.session_state['traits_df'],
            '',
            users_settings,
            {'selected_metadata': selected_metadata, 'taxonomic_level': taxonomic_level}
        )

        st.plotly_chart(fig_read_props_1, use_container_width=True, config=st.session_state['config'])
        export_plot(st.session_state['TaXon_table_xlsx'], 'Read_proportions_plots', 'heatmap', selected_metadata, fig_read_props_1, users_settings, 'plotly')