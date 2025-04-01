import pandas as pd
import numpy as np
from pathlib import Path
from taxontabletools2.utilities import filter_taxontable
import plotly.graph_objects as go
import plotly.express as px
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
import streamlit as st
import plotly.colors

def pcoa_analysis(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, selected_metadata, traits_df, selected_traits, user_settings, tool_settings):

    # Make copies of input dataframes to prevent altering the originals
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    # Use order of samples as sorted in the metadata_df
    samples = metadata_df['Sample'].values.tolist()

    # Extract relevant settings from the tool settings
    selected_metadata = tool_settings['selected_metadata']
    taxonomic_level = tool_settings['taxonomic_level']
    metric = tool_settings['metric']

    if selected_metadata != 'All samples':

        ## Select only the relevant samples
        sub_df = metadata_df[['Sample', selected_metadata]]
        sub_df = sub_df.loc[sub_df[selected_metadata] != '']
        samples = sub_df['Sample'].values.tolist()
        sample_categories = sub_df[selected_metadata].values.tolist()

        # Create colors for categories
        colorscale_name = user_settings.get('colorsequence', 'Plotly')
        colors = getattr(px.colors.qualitative, colorscale_name, px.colors.qualitative.Plotly)*200
        categories = sub_df[selected_metadata].drop_duplicates().values.tolist()
        color_dict = {i:j for i,j in zip(categories, colors)}
        sample_colors = [color_dict[metadata] for sample, metadata in sub_df.values.tolist()]

    else:
        # Create a color for each sample (use continous colors here)
        sample_categories = samples
        num_colors = len(sample_categories)
        sample_colors = plotly.colors.sample_colorscale(user_settings['colorscale'], [i / (num_colors - 1) for i in range(num_colors)])

    # Get unique taxa at the specified taxonomic level and exclude empty values
    all_taxa = [i for i in taxon_table_df[taxonomic_level].dropna().unique() if i != '']

    # Initialize result dictionary
    res = {}

    # Case where all samples are used
    for sample in samples:
        # Sum the reads for each taxon in the current sample and filter for non-zero values
        sample_taxa = taxon_table_df[[taxonomic_level, sample]].groupby(taxonomic_level).sum()
        sample_taxa = sample_taxa.loc[sample_taxa[sample] > 0, sample]

        # Store binary presence/absence data if metric is 'Jaccard', otherwise store actual values
        if metric == 'Jaccard':
            res[sample] = [1 if taxon in sample_taxa.index else 0 for taxon in all_taxa]
        else:
            res[sample] = [sample_taxa[taxon] if taxon in sample_taxa.index else 0 for taxon in all_taxa]

    # Convert result dict into DataFrame (rows: taxa, columns: samples)
    presence_absence_df = pd.DataFrame(res, index=all_taxa).transpose()

    # Calculate distance between samples
    if metric == 'Jaccard':
        distances = pdist(presence_absence_df, metric=metric)
    else:
        distances = pdist(presence_absence_df, metric=metric)

    # Convert to square form distance matrix
    distance_matrix = squareform(distances)

    # Convert to a DataFrame for easier manipulation
    distance_matrix_df = pd.DataFrame(distance_matrix, index=samples, columns=samples)

    # Try to import skbio (which is difficult to install)
    try:
        from skbio.stats.ordination import pcoa
        from skbio import DistanceMatrix
        from skbio.stats.distance import anosim
        skbio = True
    except ImportError:
        skbio = False

    # Only continue if skbio is installed
    if skbio == True:

        # Calculate pcoa
        pcoa_res = pcoa(distance_matrix_df)

        # Collect expained variance
        pcoa_explained_variance_df = pd.DataFrame(pcoa_res.proportion_explained, columns=['explained_variance'])*100
        pcoa_explained_variance_df = pcoa_explained_variance_df[pcoa_explained_variance_df['explained_variance'] > 1]

        # Collect values and already filter for PC axes with >1% explained variance
        pcoa_df = pd.DataFrame(pcoa_res.samples)
        pcoa_df.index=samples
        pcoa_axes = list(pcoa_explained_variance_df.index)
        pcoa_df = pcoa_df[pcoa_axes]
        pcoa_df['Color'] = sample_colors
        pcoa_df['Metadata'] = sample_categories

        # Create a dict for display and to draw from later
        axes_dict = pcoa_explained_variance_df.to_dict()['explained_variance']

        # Calculate anosim
        # Calculate ANOSIM if specific metadata is selected
        anosim_result = None
        if selected_metadata != 'All samples':
            # Create DistanceMatrix object for ANOSIM
            dm = DistanceMatrix(distance_matrix)

            # Perform ANOSIM using the selected metadata as grouping variable
            result = anosim(dm, sample_categories, permutations=999)
            anosim_result = result  # Store the ANOSIM result

        return pcoa_df, axes_dict, anosim_result, selected_metadata

    else:
        st.error('Please install skbio to calculate PCoAs!')

def draw_outlines(fig, x_values, y_values, metadata, color):
    ## collect samples that form the outline
    x_plane, y_plane = [], []
    hull = ConvexHull(np.column_stack((x_values, y_values)))
    for i in hull.vertices:
        x_plane.append(x_values[i])
        y_plane.append(y_values[i])
    x_plane.append(x_values[hull.vertices[0]])
    y_plane.append(y_values[hull.vertices[0]])

    ## draw the outline
    fig.add_trace(go.Scatter(x=x_plane, y=y_plane, mode='lines', name=metadata, marker_color=color, fill='toself'))

def display_pcoa(pcoa_df, axes_dict, selected_axes, user_settings, tool_settings):

        # Extract axes to display
        display_axes = [i.split(' ')[0] for i in selected_axes if i != 'Not selected']
        explained_variances = {i:axes_dict[i] for i in display_axes}
        pcoa_axes = list(explained_variances.keys())
        pcoa_df_selected_axes = pcoa_df[pcoa_axes]

        ## Display anosim
        anosim_result = tool_settings['anosim_result']
        if anosim_result is not None:
            r_value = anosim_result['test statistic']
            p_value = anosim_result['p-value']
            title = f'{tool_settings["pcoa_metadata"]}; {tool_settings["taxonomic_level"]}; ANOSIM R={r_value}, p={p_value}'
        else:
            title = f'{tool_settings["pcoa_metadata"]}; {tool_settings["taxonomic_level"]}; ANOSIM not calculated'

        if len(explained_variances) == 2:
            pcoa_2D_df = pcoa_df_selected_axes.copy()
            pcoa_2D_df['Color'] = pcoa_df['Color'].values.tolist()
            pcoa_2D_df['Metadata'] = pcoa_df['Metadata'].values.tolist()

            fig = go.Figure()

            for metadata in set(pcoa_df['Metadata'].values.tolist()):
                sub_df = pcoa_2D_df.loc[pcoa_2D_df['Metadata'] == metadata]

                x_values = sub_df[pcoa_axes[0]].values.tolist()
                y_values = sub_df[pcoa_axes[1]].values.tolist()
                text_values = list(sub_df.index)
                colors = sub_df['Color'].values.tolist()

                fig.add_trace(go.Scatter(x=x_values, y=y_values, text=text_values, marker_color=colors, marker=dict(size=user_settings['scatter_size']),  name=metadata, mode='markers'))

            if tool_settings['pcoa_metadata'] != 'All samples' and tool_settings['draw_outlines'] == True:
                for category in set(pcoa_df['Metadata'].values.tolist()):
                    sub_df = pcoa_df.loc[pcoa_df['Metadata'] == category]
                    if len(sub_df) >= 2:
                        x_values = sub_df[pcoa_axes[0]].values.tolist()
                        y_values = sub_df[pcoa_axes[1]].values.tolist()
                        color = sub_df['Color'].drop_duplicates().values.tolist()[0]
                        draw_outlines(fig, x_values, y_values, category, color)

            # Update layout
            fig.update_xaxes(title=f'{selected_axes[0]}')
            fig.update_yaxes(title=f'{selected_axes[1]}')
            fig.update_layout(title=title,
                              template=user_settings['template'],
                              font_size=user_settings['font_size'],
                              width=user_settings['plot_width'],
                              height=user_settings['plot_height'],
                              showlegend=user_settings['show_legend'],
                              yaxis_title_font=dict(size=user_settings['font_size']),
                              yaxis_tickfont=dict(size=user_settings['font_size']),
                              xaxis_title_font=dict(size=user_settings['font_size']),
                              xaxis_tickfont=dict(size=user_settings['font_size']))

            return fig

        elif len(explained_variances) == 3:

            pcoa_3D_df = pcoa_df_selected_axes.copy()
            pcoa_3D_df['Color'] = pcoa_df['Color'].values.tolist()
            pcoa_3D_df['Metadata'] = pcoa_df['Metadata'].values.tolist()

            fig = go.Figure()

            for metadata in set(pcoa_df['Metadata'].values.tolist()):
                sub_df = pcoa_3D_df.loc[pcoa_3D_df['Metadata'] == metadata]

                x_values = sub_df[pcoa_axes[0]].values.tolist()
                y_values = sub_df[pcoa_axes[1]].values.tolist()
                z_values = sub_df[pcoa_axes[2]].values.tolist()
                text_values = list(sub_df.index)
                colors = sub_df['Color'].values.tolist()

                fig.add_trace(go.Scatter3d(x=x_values, y=y_values, z=z_values, text=text_values, marker_color=colors, marker=dict(size=user_settings['scatter_size']),  name=metadata, mode='markers'))

            # Update layout with custom axis titles
            fig.update_layout(
                scene=dict(
                    xaxis_title=f'{selected_axes[0]}',  # Replace with your X axis title
                    yaxis_title=f'{selected_axes[1]}',  # Replace with your Y axis title
                    zaxis_title=f'{selected_axes[2]}'  # Replace with your Z axis title
                )
            )

            fig.update_layout(title=title,
                              template=user_settings['template'],
                              font_size=user_settings['font_size'],
                              width=user_settings['plot_width'],
                              height=user_settings['plot_height'],
                              showlegend=user_settings['show_legend'],
                              yaxis_title_font=dict(size=user_settings['font_size']),
                              yaxis_tickfont=dict(size=user_settings['font_size']),
                              xaxis_title_font=dict(size=user_settings['font_size']),
                              xaxis_tickfont=dict(size=user_settings['font_size']))

            return fig






