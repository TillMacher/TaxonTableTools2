## variables for testing: IGNORE
path_to_outdirs = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test/')
taxon_table_xlsx = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test/TaXon_tables/GeDNA_MZB_20_21_taxon_table_cons_NCsub_MZB_merged_shared.xlsx')
taxon_table_df = pd.read_excel(taxon_table_xlsx).fillna('')
taxon_table_df = strip_traits(taxon_table_df)
samples = list(taxon_table_df.columns[9:])
metadata_df = pd.read_excel(taxon_table_xlsx, sheet_name='Metadata Table').fillna('')
selected_metadata = 'Venn3'
traits_df = collect_traits(taxon_table_df)
selected_traits = ''
user_settings = {
    'project_name': 'AA_test',
    'path_to_outdirs': path_to_outdirs,
    'plot_height': 800,
    'plot_width': 800,
    'template': 'simple_white',
    'font_size': 16,
    'clustering_unit': 'ESVs',
    'scatter_size': 10,
    'color_1': 'Navy',
    'color_2': 'Teal',
    'colorsequence': 'Plotly',
    'colorscale': 'Blues',
    'show_legend':True,
    'show_xaxis':True,
    'show_yaxis':False,
    }
tool_settings = {'selected_metadata': 'Venn2', 'taxonomic_level': 'Species', 'metric': 'Jaccard', 'dimensions':'2', 'draw_outlines':True, }
# taxon_table_df = simple_taxontable(taxon_table_xlsx, taxon_table_df, samples, metadata_df, False)


selected_axes = ['PC1 (57.86%)', 'PC2 (14.89%)', 'PC3 (13.19%)']

# Use color map based on user settings
colorscale_name = user_settings.get('colorsequence', 'Plotly')
colors = getattr(px.colors.qualitative, colorscale_name, px.colors.qualitative.Plotly)

# Use color map based on user settings
colorscale_name = user_settings.get('colorscale', 'Viridis')

tool_settings['perlodes_TTT_conversion_xlsx'] = "/Library/Frameworks/Python.framework/Versions/3.10/lib/python3.10/site-packages/taxontabletools2/WFD_conversion/perlodes_TTT_conversion.xlsx"
tool_settings['selected_metadata'] = 'Type'
tool_settings['presence_absence'] = True

import plotly.colors

# Sample 100 colors from the colorscale
num_colors = 100
color_list = plotly.colors.sample_colorscale(user_settings['colorscale'], [i / (num_colors - 1) for i in range(num_colors)])
