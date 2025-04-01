import streamlit as st
import pandas as pd
import os
import sys
import subprocess
import webbrowser
import glob
from pathlib import Path
from taxontabletools2.utilities import load_df, collect_traits, strip_traits
import plotly.express as px

def TTT_variables():

    directories_to_create = ["Venn_diagrams","TaXon_tables", "Rarefaction_curves",
                             "Site_occupancy_plots", "Read_proportions_plots",
                             "Krona_charts", "Alpha_diversity", "Beta_diversity",
                             "PCoA_plots", "Replicate_analysis", "GBIF", "Occurrence_analysis",
                             "Per_taxon_statistics", "NMDS_plots", "Table_comparison", "Perlodes",
                             "EFI", "Diathor", "Phylib", "Time_series", "Basic_stats", "Fasta", "Import",
                             "Rarefaction_curves"]

    available_templates_list = ['seaborn', 'ggplot2', 'simple_white', 'plotly', 'plotly_dark', 'presentation', 'plotly_white']

    available_clustering_units = ['OTUs', 'zOTUs', 'ESVs', 'ASVs']

    plotly_colors = ["aliceblue", "antiquewhite", "aqua", "aquamarine", "azure",
    "beige", "bisque", "black", "blanchedalmond", "blue",
    "blueviolet", "brown", "burlywood", "cadetblue",
    "chartreuse", "chocolate", "coral", "cornflowerblue",
    "cornsilk", "crimson", "cyan", "darkblue", "darkcyan",
    "darkgoldenrod", "darkgray", "darkgrey", "darkgreen",
    "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange",
    "darkorchid", "darkred", "darksalmon", "darkseagreen",
    "darkslateblue", "darkslategray", "darkslategrey",
    "darkturquoise", "darkviolet", "deeppink", "deepskyblue",
    "dimgray", "dimgrey", "dodgerblue", "firebrick",
    "floralwhite", "forestgreen", "fuchsia", "gainsboro",
    "ghostwhite", "gold", "goldenrod", "gray", "grey", "green",
    "greenyellow", "honeydew", "hotpink", "indianred", "indigo",
    "ivory", "khaki", "lavender", "lavenderblush", "lawngreen",
    "lemonchiffon", "lightblue", "lightcoral", "lightcyan",
    "lightgoldenrodyellow", "lightgray", "lightgrey",
    "lightgreen", "lightpink", "lightsalmon", "lightseagreen",
    "lightskyblue", "lightslategray", "lightslategrey",
    "lightsteelblue", "lightyellow", "lime", "limegreen",
    "linen", "magenta", "maroon", "mediumaquamarine",
    "mediumblue", "mediumorchid", "mediumpurple",
    "mediumseagreen", "mediumslateblue", "mediumspringgreen",
    "mediumturquoise", "mediumvioletred", "midnightblue",
    "mintcream", "mistyrose", "moccasin", "navajowhite", "navy",
    "oldlace", "olive", "olivedrab", "orange", "orangered",
    "orchid", "palegoldenrod", "palegreen", "paleturquoise",
    "palevioletred", "papayawhip", "peachpuff", "peru", "pink",
    "plum", "powderblue", "purple", "red", "rosybrown",
    "royalblue", "saddlebrown", "salmon", "sandybrown",
    "seagreen", "seashell", "sienna", "silver", "skyblue",
    "slateblue", "slategray", "slategrey", "snow", "springgreen",
    "steelblue", "tan", "teal", "thistle", "tomato", "turquoise",
    "violet", "wheat", "white", "whitesmoke", "yellow",
    "yellowgreen"]

    available_colorsequences = {"Plotly":px.colors.qualitative.Plotly, "G10":px.colors.qualitative.G10,
    "T10":px.colors.qualitative.T10, "Alphabet":px.colors.qualitative.Alphabet, "Dark24":px.colors.qualitative.Dark24, "Dark24_r":px.colors.qualitative.Dark24_r,
    "Light24":px.colors.qualitative.Light24, "Set1":px.colors.qualitative.Set1, "Pastel1":px.colors.qualitative.Pastel,
    "Dark2":px.colors.qualitative.Dark2, "Set2":px.colors.qualitative.Set2, "Pastel2":px.colors.qualitative.Pastel2,
    "Set3":px.colors.qualitative.Set3, "Antique":px.colors.qualitative.Antique,"Bold":px.colors.qualitative.Bold,
    "Pastel":px.colors.qualitative.Pastel, "Prism":px.colors.qualitative.Prism, "Safe":px.colors.qualitative.Safe,
    "Vivid":px.colors.qualitative.Vivid}

    available_colorscales = px.colors.named_colorscales()

    available_taxonomic_levels_list= ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'ID']
    st.session_state['available_taxonomic_levels_list'] = available_taxonomic_levels_list

    return [directories_to_create, available_templates_list, available_clustering_units, plotly_colors, available_colorsequences, available_colorscales, available_taxonomic_levels_list]

def start():
    # Collect path to default files
    path_to_ttt = Path(__file__).resolve().parent
    user_preferences_xlsx = path_to_ttt.joinpath('user_preferences.xlsx')
    user_preferences_df = pd.read_excel(user_preferences_xlsx)
    path_to_projects = Path(user_preferences_df.loc[user_preferences_df['Variable'] == 'path_to_outdirs', 'Value'].values[0])

    projects_dict = {}
    for project in sorted(glob.glob(str(path_to_projects.joinpath('*')))):
        projects_dict[Path(project).stem] = Path(project)

    if projects_dict == {}:
        projects_dict['Default_project'] = str(path_to_projects.joinpath('Default_project'))

    # If this is the first load, initialize the session state values.
    if 'TaXon_table_df' not in st.session_state:
        st.session_state['TaXon_table_df'] = pd.DataFrame()

    if 'project_name' not in st.session_state:
        st.session_state['project_name'] = list(projects_dict.keys())[0]

    directories_to_create, available_templates_list, available_clustering_units, plotly_colors, available_colorsequences, available_colorscales, available_taxonomic_levels_list = TTT_variables()

    user_preferences_dict = {key: value for (key, value) in user_preferences_df.values.tolist()}
    default_template = user_preferences_dict['template']
    default_height = user_preferences_dict['plot_height']
    default_width = user_preferences_dict['plot_width']
    default_show_legend = user_preferences_dict['show_legend']
    default_color1 = user_preferences_dict['color_1']
    default_color2 = user_preferences_dict['color_2']
    default_colorsequence = user_preferences_dict['colorsequence']
    default_colorscale = user_preferences_dict['colorscale']
    default_font_size = int(user_preferences_dict['font_size'])
    default_clustering_unit = user_preferences_dict['clustering_unit']

    # Create empty df to store intermediate TaXon tables
    cache_df = pd.DataFrame()
    TaXon_table_df = pd.DataFrame()

    ### Sidebar

    st.sidebar.write(""" # Project """)
    st.sidebar.write(""" ### TaXon Table """)

    ## CALL PROJECTS
    available_projects_list = list(projects_dict.keys())
    if 'project_name' not in st.session_state:
        st.session_state['project_name'] = st.sidebar.selectbox('Select a project:', available_projects_list)
        st.session_state['path_to_outdirs'] = projects_dict[st.session_state['project_name']]
    else:
        value = st.session_state['project_name']
        st.session_state['project_name'] = st.sidebar.selectbox('Select a project:', available_projects_list, index=available_projects_list.index(value))
        st.session_state['path_to_outdirs'] = projects_dict[st.session_state['project_name']]

    ## CALL PROJECTS
    if st.sidebar.button('Reveal project in finder'):
        if sys.platform == 'win32':
            subprocess.Popen(['start', st.session_state['path_to_outdirs']], shell=True)
        elif sys.platform == 'darwin':
            subprocess.Popen(['open', st.session_state['path_to_outdirs']])

    ## Make sure all directories exist
    [os.makedirs(Path(st.session_state['path_to_outdirs']).joinpath(i), exist_ok=True) for i in directories_to_create]

    ## load TaXon table
    TaXon_table_upload = st.sidebar.file_uploader('Load TaXon table')

    if 'TaXon_table_df' not in st.session_state:
        st.session_state['TaXon_table_df'] = pd.DataFrame()

    if st.sidebar.button('Load table'):
        if TaXon_table_upload is not None:
            # Reset the session state before loading new data
            st.session_state['TaXon_table_df'] = pd.DataFrame()  # Clear the old table

            # Load the new table
            TaXon_table_xlsx = Path('{}/TaXon_tables/{}'.format(st.session_state['path_to_outdirs'], TaXon_table_upload.name))
            st.session_state['TaXon_table_xlsx'] = TaXon_table_xlsx

            # Verify that the file exists before attempting to load it
            if not st.session_state['TaXon_table_xlsx'].is_file():
                error1, error2, error3 = st.columns(3)
                with error2:
                    st.write('##### The TaXon table does not exist in the selected project folder!')
            else:
                # Load the TaXon table
                TaXon_table_df = load_df(st.session_state['TaXon_table_xlsx'])

                # Update session state with the new table
                st.session_state['TaXon_table_df'] = TaXon_table_df

                # Process the traits and other data
                traits_df = collect_traits(st.session_state['TaXon_table_df'])
                st.session_state['traits_df'] = traits_df
                TaXon_table_df = strip_traits(st.session_state['TaXon_table_df'])
                st.session_state['TaXon_table_df'] = TaXon_table_df
                samples = list(st.session_state['TaXon_table_df'].columns[9:])
                st.session_state['samples'] = samples

                # Load the metadata if it exists
                try:
                    metadata_df = pd.read_excel(st.session_state['TaXon_table_xlsx'], sheet_name='Metadata Table').fillna('')
                    st.session_state['metadata_df'] = metadata_df
                    if sorted(st.session_state['metadata_df']['Sample'].values.tolist()) != sorted(st.session_state['samples']):
                        st.session_state['metadata_df'] = pd.DataFrame(columns=['Sample', 'Metadata'])
                except ValueError:
                    dummy_metadata = [[j,i] for i,j in enumerate(samples)]
                    st.session_state['metadata_df'] = pd.DataFrame(dummy_metadata, columns=['Sample', 'Placeholder'])

                if sorted(samples) != sorted(st.session_state['metadata_df']['Sample'].values.tolist()):
                    dummy_metadata = [[j,i] for i,j in enumerate(samples)]
                    st.session_state['metadata_df'] = pd.DataFrame(dummy_metadata, columns=['Sample', 'Placeholder'])
                    st.error('Warning: The samples in the metadata sheet do not match with the samples in the taxontable sheet!')


        else:
            st.warning('Please provide a Taxon Table file first!')

    if 'TaXon_table_xlsx' in st.session_state:
        st.sidebar.write(Path(st.session_state['TaXon_table_xlsx']).name)

        if st.sidebar.button('Open table'):
            file_path = str(Path(st.session_state['TaXon_table_xlsx']))
            if sys.platform.startswith("win"):  # Windows
                os.startfile(file_path)
            elif sys.platform.startswith("darwin"):  # macOS
                subprocess.run(["open", file_path])
            else:  # Linux
                subprocess.run(["xdg-open", file_path])

    ####################################################################################################################

    st.sidebar.write(""" # Settings """)
    st.sidebar.write(""" ### Presentation """)

    # PLOT HEIGHT
    if 'plot_height' not in st.session_state:
        st.session_state['plot_height'] = st.sidebar.number_input('Height', 400, 4000, default_height)
    else:
        value = st.session_state['plot_height']
        st.session_state['plot_height'] = st.sidebar.number_input('Height', 400, 4000, value)

    # PLOT WIDTH
    if 'plot_width' not in st.session_state:
        st.session_state['plot_width'] = st.sidebar.number_input('Width', 400, 4000, default_width)
    else:
        value = st.session_state['plot_width']
        st.session_state['plot_width'] = st.sidebar.number_input('Width', 400, 4000, value)

    # SHOW LEGEND
    if 'show_legend' not in st.session_state:
        st.session_state['show_legend'] = st.sidebar.selectbox('Show legend', [True, False], index=[True, False].index(default_show_legend))
    else:
        value = st.session_state['show_legend']
        st.session_state['show_legend'] = st.sidebar.selectbox('Show legend', [True, False], index=[True, False].index(value))

    # TEMPLATE
    if 'template' not in st.session_state:
        st.session_state['template'] = st.sidebar.selectbox('Layout', available_templates_list, index=available_templates_list.index(default_template))
    else:
        value = st.session_state['template']
        st.session_state['template'] = st.sidebar.selectbox('Layout', available_templates_list,
                                                            index=available_templates_list.index(value))

    # FONT_SIZE
    if 'font_size' not in st.session_state:
        st.session_state['font_size'] = st.sidebar.slider('Font size', 6, 30, default_font_size)
    else:
        value = st.session_state['font_size']
        st.session_state['font_size'] = st.sidebar.slider('Font size', 6, 30, value)

    # CLUSTERING UNIT
    if 'clustering_unit' not in st.session_state:
        st.session_state['clustering_unit'] = st.sidebar.selectbox('Clustering Unit', available_clustering_units,
                                                                   index=available_clustering_units.index(
                                                                       default_clustering_unit))
    else:
        value = st.session_state['clustering_unit']
        st.session_state['clustering_unit'] = st.sidebar.selectbox('Clustering Unit', available_clustering_units,
                                                                   index=available_clustering_units.index(value))

    # SCATTER SIZE
    if 'scatter_size' not in st.session_state:
        st.session_state['scatter_size'] = st.sidebar.number_input('Scatter size', 0, 40, 15)
    else:
        value = st.session_state['scatter_size']
        st.session_state['scatter_size'] = st.sidebar.number_input('Scatter size', 0, 40, value)

    st.sidebar.write(""" ### Colors """)

    # COLOR 1
    if 'color_1' not in st.session_state:
        st.session_state['color_1'] = st.sidebar.selectbox('Color 1', plotly_colors,
                                                           index=plotly_colors.index(default_color1))
    else:
        value = st.session_state['color_1']
        st.session_state['color_1'] = st.sidebar.selectbox('Color 1', plotly_colors, index=plotly_colors.index(value))

    # COLOR 2
    if 'color_2' not in st.session_state:
        st.session_state['color_2'] = st.sidebar.selectbox('Color 2', plotly_colors,
                                                           index=plotly_colors.index(default_color2))
    else:
        value = st.session_state['color_2']
        st.session_state['color_2'] = st.sidebar.selectbox('Color 2', plotly_colors, index=plotly_colors.index(value))

    # COLORSEQUENCE
    if 'colorsequence' not in st.session_state:
        st.session_state['colorsequence'] = st.sidebar.selectbox('Color sequence',
                                                                 list(available_colorsequences.keys()),
                                                                 index=list(available_colorsequences.keys()).index(
                                                                     default_colorsequence))
    else:
        value = st.session_state['colorsequence']
        st.session_state['colorsequence'] = st.sidebar.selectbox('Color sequence',
                                                                 list(available_colorsequences.keys()),
                                                                 index=list(available_colorsequences.keys()).index(
                                                                     value))

    # COLORSCALE
    if 'colorscale' not in st.session_state:
        st.session_state['colorscale'] = st.sidebar.selectbox('Color scale', available_colorscales,
                                                              index=available_colorscales.index(default_colorscale))
    else:
        value = st.session_state['colorscale']
        st.session_state['colorscale'] = st.sidebar.selectbox('Color scale', available_colorscales,
                                                              index=available_colorscales.index(value))

    new_user_preferences_df = pd.DataFrame([['project_name', st.session_state['project_name']],
                                            ['path_to_outdirs', st.session_state['path_to_outdirs']],
                                            ['plot_height', st.session_state['plot_height']],
                                            ['plot_width', st.session_state['plot_width']],
                                            ['template', st.session_state['template']],
                                            ['show_legend', st.session_state['show_legend']],
                                            ['font_size', st.session_state['font_size']],
                                            ['clustering_unit', st.session_state['clustering_unit']],
                                            ['scatter_size', st.session_state['scatter_size']],
                                            ['color_1', st.session_state['color_1']],
                                            ['color_2', st.session_state['color_2']],
                                            ['colorsequence', st.session_state['colorsequence']],
                                            ['colorscale', st.session_state['colorscale']],
                                            ], columns=['Variable', 'Value'])

    new_user_preferences_dict = {
        'project_name': st.session_state['project_name'],
        'path_to_outdirs': path_to_projects,
        'plot_height': st.session_state['plot_height'],
        'plot_width': st.session_state['plot_width'],
        'show_legend': st.session_state['show_legend'],
        'template': st.session_state['template'],
        'font_size': st.session_state['font_size'],
        'clustering_unit': st.session_state['clustering_unit'],
        'scatter_size': st.session_state['scatter_size'],
        'color_1': st.session_state['color_1'],
        'color_2': st.session_state['color_2'],
        'colorsequence': st.session_state['colorsequence'],
        'colorscale': st.session_state['colorscale']
    }

    st.sidebar.write(""" ### Save """)

    if st.sidebar.button('Save preferences'):
        updated_user_preferences_df = pd.DataFrame(list(new_user_preferences_dict.items()), columns=['Variable', 'Value'])
        updated_user_preferences_df.to_excel(user_preferences_xlsx, index=False)
        st.sidebar.write('Saved preferences...!')

    if st.sidebar.button('Help'):
        webbrowser.open_new_tab(
            'https://sites.google.com/d/1JHe4k_-j7X0_XHj0eQ9HRJYfLwb-sVt4/p/1mSo5PZ0O5oYSz6jLGLHS5wAKyw0ZBo_C/edit')

    return {
        'projects_dict': projects_dict,
        'cache_df': cache_df,
        'TaXon_table_df': TaXon_table_df,
        'new_user_preferences_dict': new_user_preferences_dict
        }
