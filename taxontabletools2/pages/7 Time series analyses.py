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
from taxontabletools2.utilities import export_plot
from taxontabletools2.table_processing import replicate_merging
from taxontabletools2.table_processing import negative_control_subtraction
from taxontabletools2.table_processing import read_based_filter
from taxontabletools2.table_processing import read_based_normalisation
from taxontabletools2.table_processing import taxonomic_filtering
from taxontabletools2.table_processing import sample_filtering
from taxontabletools2.time_series import time_series_richness_with_ci
import streamlit as st
from taxontabletools2.start import start

import numpy as np
import pymannkendall as mk

# Data generation for analysis
data = np.random.rand(360,1)

result = mk.original_test(data)
print(result)

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## Time series """)

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
    st.write('### Richness')

    time1, time2, time3, time4, time5 = st.columns(5)

    with time1:
        available_metadata = ['All samples'] + st.session_state['metadata_df'].columns.tolist()[1:]
        selected_metadata = st.selectbox('Select metadata to test:', available_metadata)

    with time2:
        taxonomic_level_1 = st.selectbox('Select first taxonomic level:', st.session_state['available_taxonomic_levels_list'], index=5)

    with time3:
        taxonomic_level_2 = st.selectbox('Select second taxonomic level:', st.session_state['available_taxonomic_levels_list'], index=3)

    with time4:
        n_bootstraps = st.number_input('Select number of LOESS bootstraps:', min_value=0, max_value=10000, value=1000)

    with time5:
        ci_level = st.number_input('Select confidence intervall :', min_value=0, max_value=100, value=95)

    if st.button(f'Calculate {taxonomic_level_1} richness diagram for {selected_metadata}'):

        # Call the venn function and render it within Streamlit
        time_series_richness_with_ci_fig = time_series_richness_with_ci(
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
             'taxonomic_level_1': taxonomic_level_1,
             'taxonomic_level_2': taxonomic_level_2,
             'n_bootstraps': n_bootstraps,
             'ci_level': ci_level
             }
        )

        st.plotly_chart(time_series_richness_with_ci_fig, use_container_width=True, config=st.session_state['config'])
        export_plot(st.session_state['TaXon_table_xlsx'], 'Time_series', 'richness', selected_metadata, time_series_richness_with_ci_fig, users_settings, 'plotly')