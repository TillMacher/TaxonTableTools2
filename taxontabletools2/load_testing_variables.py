from taxontabletools2.utilities import load_df, collect_traits, strip_traits
from pathlib import Path
import pandas as pd
import os

# DiATHOR TEST
path_to_outdirs = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test')
taxon_table_xlsx = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test/TaXon_tables/LL_diatoms_taxon_table_cons_NCsub_merged_diatoms.xlsx')
taxon_table_df = load_df(taxon_table_xlsx)
selected_metadata = 'River type'
traits_df = collect_traits(taxon_table_df)
taxon_table_df = strip_traits(taxon_table_df)
samples = list(taxon_table_df.columns[9:])
metadata_df = pd.read_excel(taxon_table_xlsx, sheet_name='Metadata Table').fillna('')
user_settings = {'project_name': 'AA_test', 'path_to_outdirs': Path('/Users/tillmacher/Desktop/TTT_projects'),
                 'plot_height': 1000, 'plot_width': 1000, 'show_legend': True, 'template': 'simple_white',
                 'font_size': 20, 'clustering_unit': 'ESVs', 'scatter_size': 15, 'color_1': 'navy', 'color_2': 'teal',
                 'colorsequence': 'Plotly', 'colorscale': 'blues'}
users_settings = user_settings
tool_settings = {'selected_metadata': 'River type', 'presence_absence': True}
perlodes_TTT_conversion_xlsx = Path('/Users/tillmacher/Documents/GitHub/TaxonTableTools2.0/taxontabletools2/WFD_conversion/perlodes_TTT_conversion.xlsx')
script_dir = "/Users/tillmacher/Applications/miniconda3/miniconda3/envs/TTT/lib/python3.12/site-packages/taxontabletools2"
script_path = Path(os.path.join(script_dir, 'Rscripts',  'Diathor.R'))





# EFI TEST
path_to_outdirs = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test')
taxon_table_xlsx = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test/TaXon_tables/LL_vertebrates_taxon_table_cons_NCsub_merged_vertebrates_fish.xlsx')
taxon_table_df = load_df(taxon_table_xlsx)
traits_df = collect_traits(taxon_table_df)
taxon_table_df = strip_traits(taxon_table_df)
samples = list(taxon_table_df.columns[9:])
metadata_df = pd.read_excel(taxon_table_xlsx, sheet_name='Metadata Table').fillna('')
user_settings = {'project_name': 'AA_test', 'path_to_outdirs': Path('/Users/tillmacher/Desktop/TTT_projects'),
                 'plot_height': 1000, 'plot_width': 1000, 'show_legend': True, 'template': 'simple_white',
                 'font_size': 20, 'clustering_unit': 'ESVs', 'scatter_size': 15, 'color_1': 'navy', 'color_2': 'teal',
                 'colorsequence': 'Plotly', 'colorscale': 'blues'}
users_settings = user_settings
efi_columns = [
    "EFI_Day", "EFI_Month", "EFI_Year", "EFI_Longitude", "EFI_Latitude",
    "EFI_Actual.river.slope", "EFI_Temp.jul", "EFI_Temp.jan", "EFI_Floodplain.site",
    "EFI_Water.source.type", "EFI_Geomorph.river.type", "EFI_Distance.from.source",
    "EFI_Area.ctch", "EFI_Natural.sediment", "EFI_Ecoreg", "EFI_Eft.type",
    "EFI_Fished.area", "EFI_Method"
]
tool_settings = {'efi_columns': efi_columns, 'presence_absence': True}
script_dir = "/Users/tillmacher/Applications/miniconda3/miniconda3/envs/TTT/lib/python3.12/site-packages/taxontabletools2"
script_path = Path(os.path.join(script_dir, 'Rscripts',  'EFI.R'))
higherTaxonKey = 1


# Rarefaction
tool_settings = {'reps': 1000, 'taxonomic_level': 'Species'}


# PERLODES TEST
path_to_outdirs = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test')
taxon_table_xlsx = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test/TaXon_tables/LL_invertebrates_taxon_table_cons_NCsub_merged_invertebrates_aquatic.xlsx')
taxon_table_df = load_df(taxon_table_xlsx)
selected_metadata = 'River type'
traits_df = collect_traits(taxon_table_df)
taxon_table_df = strip_traits(taxon_table_df)
samples = list(taxon_table_df.columns[9:])
metadata_df = pd.read_excel(taxon_table_xlsx, sheet_name='Metadata Table').fillna('')
user_settings = {'project_name': 'AA_test', 'path_to_outdirs': Path('/Users/tillmacher/Desktop/TTT_projects'),
                 'plot_height': 1000, 'plot_width': 1000, 'show_legend': True, 'template': 'simple_white',
                 'font_size': 20, 'clustering_unit': 'ESVs', 'scatter_size': 15, 'color_1': 'navy', 'color_2': 'teal',
                 'colorsequence': 'Plotly', 'colorscale': 'blues'}
users_settings = user_settings
tool_settings = {'selected_metadata': 'River type', 'presence_absence': True}




















