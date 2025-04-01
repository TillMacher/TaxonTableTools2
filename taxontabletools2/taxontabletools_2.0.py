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
import streamlit as st
from taxontabletools2.start import start

st.session_state['config'] = {
    'toImageButtonOptions': {'format': 'svg', 'filename': 'custom_TTT_image'},
    'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian'],
    'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'drawclosedpath', 'drawcircle', 'drawrect', 'eraseshape']
}

st.set_page_config(layout="wide")

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

# define the script location
st.session_state['script_dir'] = os.path.dirname(os.path.abspath(__file__))

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

else:
    # Display the taxon table if loaded
    st.write("""# TaxonTableTools 2.0""")

    if not st.session_state['TaXon_table_df'].empty:
        if st.session_state['TaXon_table_xlsx'].is_file() != False:
            st.write(""" ## Basic statistics """)

            ## collect stats
            st.session_state['stats_df'] = collect_sample_stats(st.session_state['TaXon_table_df'],
                                                                st.session_state['TaXon_table_xlsx'],
                                                                st.session_state['traits_df'],
                                                                st.session_state['samples'],
                                                                st.session_state['metadata_df'])

            ## PLOTS
            a1, a2 = st.columns(2)
            with a1:
                ## calculate read numbers
                fig_a1, min_reads, max_reads, avg_reads, stdev_reads = basic_stats_reads(users_settings, st.session_state['TaXon_table_df'])
                ## plot reads
                st.plotly_chart(fig_a1, use_container_width=True, config=st.session_state['config'])
                st.write('The samples contain between **{}** and **{}** reads with an average of **{} (+-{})**.'.format(
                    min_reads, max_reads, avg_reads, stdev_reads))

            with a2:
                ## calculate number of OTUs
                fig_a2, min_OTUs, max_OTUs, avg_OTUs = basic_stats_OTUs(users_settings, st.session_state['TaXon_table_df'])
                ## plot number of OTUs
                st.plotly_chart(fig_a2, use_container_width=True, config=st.session_state['config'])
                st.write('The samples contain between **{}** and **{}** {} with an average of **{}**.'.format(min_OTUs,
                                                                                                              max_OTUs,
                                                                                                              users_settings['clustering_unit'],
                                                                                                              avg_OTUs))

            a3, a4 = st.columns(2)
            with a3:
                ## calculate read numbers
                fig_a3 = taxonomic_richness(users_settings, st.session_state['TaXon_table_df'])
                ## plot reads
                st.plotly_chart(fig_a3, use_container_width=True, config=st.session_state['config'])
                st.write('Taxonomic richness')

            with a4:
                ## calculate number of OTUs
                fig_a4 = taxonomic_resolution(users_settings, st.session_state['TaXon_table_df'])
                ## plot number of OTUs
                st.plotly_chart(fig_a4, use_container_width=True, config=st.session_state['config'])
                st.write('Taxonomic resolution')

            ## Per sample information table and Metadata table
            a5, a6 = st.columns(2)
            with a5:
                st.write(""" ## Per sample information """)
                st.table(st.session_state['stats_df'])
            with a6:
                st.write(""" ## Available metadata """)
                st.table(st.session_state['metadata_df'])

        else:
            error1, error2, error3 = st.columns(3)
            with error2:
                st.write(
                    '##### Could not find your TaXon table in the selected project. Please make sure to select the according project!')

    else:
        st.error('TaXon table is empty. Please load a table to continue.')




