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
import streamlit as st
from taxontabletools2.start import start
from taxontabletools2.sample_comparison import rarefaction_curve

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## Sample comparison """)

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
    # Venn diagram
    st.write('### Venn diagram')

    venn1, venn2 = st.columns(2)

    with venn1:
        available_metadata = st.session_state['metadata_df'].columns.tolist()[1:]
        selected_metadata = st.selectbox('Select metadata to test:', available_metadata)
    with venn2:
        taxonomic_level = st.selectbox('Select a taxonomic level:', st.session_state['available_taxonomic_levels_list'], index=5)

    if st.button(f'Calculate {taxonomic_level} venn diagram for {selected_metadata}'):

        # Call the venn function and render it within Streamlit
        plt = venn(
            st.session_state['path_to_outdirs'],
            st.session_state['TaXon_table_xlsx'],
            st.session_state['TaXon_table_df'],
            st.session_state['samples'],
            st.session_state['metadata_df'],
            selected_metadata,
            st.session_state['traits_df'],
            '',
            users_settings,
            {'selected_metadata': selected_metadata,
             'taxonomic_level': taxonomic_level}
        )

        if plt:
            # Display the plot in Streamlit
            st.pyplot(plt.gcf())
            export_plot(st.session_state['TaXon_table_xlsx'], 'Venn_diagrams', 'venn', selected_metadata, plt, users_settings, 'matplot')
            plt.close()

    # Upset Diagram
    # coming soon

    ### Rarefaction curves
    st.write('Rarefaction curve')

    raf1, raf2 = st.columns(2)

    with raf1:
        repetitions_1 = st.number_input('Enter number of repetitions:', 1000, key='rarefaction_repetitions')
    with raf2:
        taxonomic_level = st.selectbox('Select a taxonomic level:', st.session_state['available_taxonomic_levels_list'], index=5, key='rarefaction_taxon')

    if st.button(f'Calculate {taxonomic_level} rarefaction curve for {selected_metadata}'):

        # Call the venn function and render it within Streamlit
        rarefaction_fig, rarefaction_df = rarefaction_curve(
            st.session_state['path_to_outdirs'],
            st.session_state['TaXon_table_xlsx'],
            st.session_state['TaXon_table_df'],
            st.session_state['samples'],
            st.session_state['metadata_df'],
            selected_metadata,
            st.session_state['traits_df'],
            '',
            users_settings,
            {'reps': repetitions_1,
             'taxonomic_level': taxonomic_level}
        )

        if rarefaction_fig:
            # Display the plot in Streamlit
            st.plotly_chart(rarefaction_fig, use_container_width=True, config=st.session_state['config'])
            export_plot(st.session_state['TaXon_table_xlsx'], 'Rarefaction_curves', 'taxon_rarefaction', selected_metadata, rarefaction_fig, users_settings, 'plotly')
