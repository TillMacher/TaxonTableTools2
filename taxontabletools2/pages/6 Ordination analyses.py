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
from taxontabletools2.diversity_analyses import alpha_boxplot
from taxontabletools2.diversity_analyses import distance_matrix
from taxontabletools2.ordination_analyses import pcoa_analysis
from taxontabletools2.ordination_analyses import display_pcoa
import streamlit as st
from taxontabletools2.start import start

# Call the sidebar function
settings = start()
user_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## Ordination analyses """)

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
    st.write('### Ordination analyses')

    st.write('### Principle Coordinate Analysis (PCoA)')

    ordination1, ordination2, ordination3 = st.columns(3)

    with ordination1:
        available_metadata = ['All samples']  + st.session_state['metadata_df'].columns.tolist()[1:]
        selected_metadata = st.selectbox('Select metadata to test:', available_metadata)

    with ordination2:
        taxonomic_level = st.selectbox('Select taxonomic level:', st.session_state['available_taxonomic_levels_list'], index=5)

    with ordination3:
        distance_metric = st.selectbox('Select metric: ', ['Jaccard', 'Braycurtis'], index=0)

    if st.button(f'Calculate {taxonomic_level} PCoA for {selected_metadata}, based on {distance_metric} similarity.'):

        st.session_state['pcoa_df'], st.session_state['axes_dict'], st.session_state['anosim_result'], st.session_state['pcoa_metadata']  = pcoa_analysis(
            st.session_state['path_to_outdirs'],
            st.session_state['TaXon_table_xlsx'],
            st.session_state['TaXon_table_df'],
            st.session_state['samples'],
            st.session_state['metadata_df'],
            selected_metadata,
            st.session_state['traits_df'],
            '',
            user_settings,
            {'selected_metadata': selected_metadata, 'taxonomic_level': taxonomic_level, 'metric': distance_metric}
        )

    if 'pcoa_df' not in st.session_state:
        st.error('Please first calculate the PCoA.')

    else:
        axis_select_1, axis_select_2, axis_select_3, axis_select_4 = st.columns(4)

        select_dict = {key:f'{key} ({round(values, 2)}%)' for key, values in st.session_state['axes_dict'].items()}

        with axis_select_1:
            axis_1 = st.selectbox('Axis 1', select_dict.values(), index=0)
        with axis_select_2:
            axis_2 = st.selectbox('Axis 2', select_dict.values(), index=1)
        with axis_select_3:
            axis_3 = st.selectbox('Axis 3', ['Not selected'] + list(select_dict.values()), index=0)
        with axis_select_4:
            draw_outlines = st.selectbox('Draw outlines according to metadata', [True, False], index=0)

        selected_axes = [axis_1, axis_2, axis_3]

        if st.button('Display PCoA'):

            fig_pcoa = display_pcoa(
                st.session_state['pcoa_df'],
                st.session_state['axes_dict'],
                selected_axes,
                user_settings,
                {'selected_metadata': selected_metadata, 'taxonomic_level': taxonomic_level, 'metric': distance_metric,
                 'draw_outlines': draw_outlines, 'anosim_result':st.session_state['anosim_result'], 'pcoa_metadata':st.session_state['pcoa_metadata']
                 }
                )

            st.plotly_chart(fig_pcoa, use_container_width=True, config=st.session_state['config'])
            export_plot(st.session_state['TaXon_table_xlsx'], 'PCoA_plots', 'PCoA', selected_metadata, fig_pcoa, user_settings, 'plotly')