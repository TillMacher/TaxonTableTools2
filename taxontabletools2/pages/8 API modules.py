import streamlit as st
import pandas as pd
import numpy as np
import subprocess
import glob, sys, time, statistics, os.path, random, datetime
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
from taxontabletools2.api_modules import gbif_accession
import streamlit as st
from taxontabletools2.start import start

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## API modules """)

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
    st.write('### GBIF | Global Biodiversity Information Facility')

    with st.expander("See more"):

        gbif1, gbif2 = st.columns(2)
        taxids_dict = {'Animalia':'1',
                          'Archae':'2',
                          'Bacteria':'3',
                          'Chromista':'4',
                          'Fungi':'5',
                          'Plantae':'6',
                          'Protozoa':'7',
                          'Viruses':'8',
                          'CUSTOM':'2436436'
                          }
        with gbif1:
            selected_higher_taxon = st.selectbox('Select a higher taxon to filter:', ['Animalia',
                                                              'Archae',
                                                              'Bacteria',
                                                              'Chromista',
                                                              'Fungi',
                                                              'Plantae',
                                                              'Protozoa',
                                                              'Viruses',
                                                              'CUSTOM'
                                                              ])
            higherTaxonKey = taxids_dict[selected_higher_taxon]

        with gbif2:
            higherTaxonKey = st.number_input('GBIF higherTaxonKey:', value=int(higherTaxonKey))

        st.write('Collected information: species key, genus key, GBIF species, synonym, IUCN red list category, and habitat')
        st.write('Note: The collected traits are added to the active TaXon table.')

        if st.button('Fetch information from GBIF'):
            gbif_accession(st.session_state['TaXon_table_xlsx'], st.session_state['TaXon_table_df'], st.session_state['samples'], st.session_state['metadata_df'], st.session_state['traits_df'], higherTaxonKey)
            current_year = datetime.datetime.now().year
            current_date = datetime.datetime.now().strftime("%d %B %Y")
            st.success('Added accession numbers to species! Please remember to reload the TaXon table!')
            st.success(f"Citation: GBIF.org ({current_year}), GBIF Home Page. Available from: https://www.gbif.org [{current_date}]")





