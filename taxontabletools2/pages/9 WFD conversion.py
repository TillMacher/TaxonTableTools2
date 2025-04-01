import webbrowser
import streamlit as st
import pandas as pd
import numpy as np
import subprocess
import glob, sys, time, statistics, os.path, random, os
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
from taxontabletools2.utilities import update_taxon_table
from taxontabletools2.table_processing import replicate_merging
from taxontabletools2.table_processing import negative_control_subtraction
from taxontabletools2.table_processing import read_based_filter
from taxontabletools2.table_processing import read_based_normalisation
from taxontabletools2.table_processing import taxonomic_filtering
from taxontabletools2.table_processing import sample_filtering
from taxontabletools2.WFD_tools import convert_to_perlodes
from taxontabletools2.WFD_tools import convert_to_phylib
from taxontabletools2.WFD_tools import convert_to_diathor
from taxontabletools2.WFD_tools import convert_to_efi
import streamlit as st
from taxontabletools2.start import start

# Call the sidebar function
settings = start()
users_settings = settings['new_user_preferences_dict']

## Page title
st.write(""" ## WFD conversion """)
st.markdown("---")

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

    ####################################################################################################################
    st.write('### 1) Macroinvertebrates')
    st.write('### Perlodes (Water Framework Directive Germany)')
    perlodes1, perlodes2 = st.columns(2)
    with perlodes1:
        presence_absence_perlodes = st.selectbox('Convert to presence/absence data:', [True, False, 'Relative', 'Genetic Diversity'], index=0, key='perlodes2')
    with perlodes2:
        if st.button('Visit the Perlodes home page to calculate status classes'):
            webbrowser.open('https://gewaesser-bewertung-berechnung.de/index.php/login.html')

    available_metadata_perlodes = st.session_state['metadata_df'].columns.tolist()[1:]
    perlodes_columns = ["Perlodes_river_type", 'Perlodes_taxa_list', 'Perlodes_usage']

    if set(perlodes_columns).issubset(set(available_metadata_perlodes)):
        st.success('All required columns to create Perlodes input file were found!')
        if st.button('Convert to Perlodes input'):
            convert_to_perlodes(st.session_state['path_to_outdirs'],
                                st.session_state['TaXon_table_xlsx'],
                                st.session_state['TaXon_table_df'],
                                st.session_state['samples'],
                                st.session_state['metadata_df'],
                                st.session_state['traits_df'],
                                {'presence_absence': presence_absence_perlodes}
                                )
    else:
        st.warning('Not all required columns are present in the metadata!')
        if st.button('Add Perlodes columns to to table'):
            for col in perlodes_columns:
                st.session_state['metadata_df'][col] = ''
            update_taxon_table(
                st.session_state['TaXon_table_xlsx'],
                st.session_state['TaXon_table_df'],
                st.session_state['traits_df'],
                st.session_state['metadata_df'],
                '')
            st.success('Required columns were added to table. Please fill out the respective Perlodes columns and reload the table!')

    with st.expander("Perlodes variables"):
        st.write('Coming soon')

    ####################################################################################################################
    st.markdown("---")
    st.write('### 2) Diatoms')
    st.write('### Phylib (Water Framework Directive Germany)')
    phylib1, phylib2 = st.columns(2)
    with phylib1:
        presence_absence_phylib = st.selectbox('Convert to presence/absence data:', [True, False, 'Relative', 'Genetic Diversity'], index=0, key='phylib2')
    with phylib2:
        if st.button('Visit the Phylib home page to calculate status classes'):
            webbrowser.open('https://gewaesser-bewertung-berechnung.de/index.php/login.html')
    if st.button(f'Create phylib upload file.'):
        print('Coming soon!')

    with st.expander("Phylib variables"):
        st.write('Coming soon')

    st.write('### Diathor (Multiple diatom indices)')
    diathor1, diathor2 = st.columns(2)
    with diathor1:
        presence_absence_diathor = st.selectbox('Convert to presence/absence data:', [True, False, 'Relative', 'Genetic Diversity'], index=0, key='diathor2')
    with diathor2:
        pass
    if st.button(f'Calculate diatom indices with Diathor.'):
        convert_to_diathor(st.session_state['path_to_outdirs'],
                           st.session_state['TaXon_table_xlsx'],
                           st.session_state['TaXon_table_df'],
                           st.session_state['samples'],
                           st.session_state['metadata_df'],
                           st.session_state['traits_df'],
                           {'presence_absence':presence_absence_diathor}
                           )

    with st.expander("Diathor variables"):
        st.write('Coming soon')

    ####################################################################################################################
    st.markdown("---")
    st.write('### 3) Fish')
    st.write('### EFI (European Fish Index)')
    efi1, efi2 = st.columns(2)
    with efi1:
        presence_absence_efi = st.selectbox('Convert to presence/absence data:', [True, False, 'Relative', 'Genetic Diversity'], index=0, key='efi2')
    with efi2:
        pass

    available_metadata_efi = st.session_state['metadata_df'].columns.tolist()[1:]
    efi_columns = [
                        "EFI_Day", "EFI_Month", "EFI_Year", "EFI_Longitude", "EFI_Latitude",
                        "EFI_Actual.river.slope", "EFI_Temp.jul", "EFI_Temp.jan", "EFI_Floodplain.site",
                        "EFI_Water.source.type", "EFI_Geomorph.river.type", "EFI_Distance.from.source",
                        "EFI_Area.ctch", "EFI_Natural.sediment", "EFI_Ecoreg", "EFI_Eft.type",
                        "EFI_Fished.area", "EFI_Method"
                    ]
    if set(efi_columns).issubset(set(available_metadata_efi)):
        st.success('All required columns to calculate the EFI were found!')
        if st.button('Calculate EFI'):
            convert_to_efi(st.session_state['path_to_outdirs'],
                               st.session_state['TaXon_table_xlsx'],
                               st.session_state['TaXon_table_df'],
                               st.session_state['samples'],
                               st.session_state['metadata_df'],
                               st.session_state['traits_df'],
                               {'efi_columns': efi_columns, 'presence_absence': presence_absence_efi}
                               )

    else:
        st.warning('Not all required columns are present in the metadata!')
        if st.button('Add EFI columns to to table'):
            for col in efi_columns:
                st.session_state['metadata_df'][col] = ''
            update_taxon_table(
                st.session_state['TaXon_table_xlsx'],
                st.session_state['TaXon_table_df'],
                st.session_state['traits_df'],
                st.session_state['metadata_df'],
                '')
            st.success('Required columns were added to table. Please fill out the respective EFI columns and reload the table!')

    with st.expander("EFI variables"):
        efi_table = Path(os.path.join(st.session_state['script_dir'], 'WFD_conversion', 'efi_table.xlsx'))
        st.table(pd.read_excel(efi_table, sheet_name='Sheet1').fillna(''))
        st.table(pd.read_excel(efi_table, sheet_name='Sheet2').fillna(''))
        efi_map = Path(os.path.join(st.session_state['script_dir'], 'WFD_conversion', 'efi_ecoregions.png'))
        caption = 'Map of Illies ecoregions represented in EFI+ and the additional Mediterranean region. Source: EFI manual'
        st.image(str(efi_map), caption=caption, use_column_width=True)