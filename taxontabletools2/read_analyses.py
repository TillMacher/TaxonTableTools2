import pandas as pd
import numpy as np
from pathlib import Path
from taxontabletools2.utilities import filter_taxontable
import plotly.graph_objects as go
import plotly.express as px

def read_props_bar(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, selected_metadata,
                   traits_df, selected_traits, user_settings, tool_settings):

    # Make copies of input dataframes to prevent altering the originals
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    # Extract relevant settings from the tool settings
    selected_metadata = tool_settings['selected_metadata']
    taxonomic_level = tool_settings['taxonomic_level']

    # Get unique taxa at the specified taxonomic level and exclude empty values
    all_taxa = taxon_table_df[taxonomic_level].dropna().unique()

    # Initialize result dictionary
    res = {}

    if selected_metadata == 'All samples':
        # Case where all samples are selected (simpler computation)
        total_reads_per_sample = taxon_table_df[samples].sum(axis=0)

        for taxon in all_taxa:
            # Filter taxon-specific data only once per taxon
            sub_df = taxon_table_df[taxon_table_df[taxonomic_level] == taxon]

            # Sum reads per taxon for each sample and calculate the proportion
            taxon_reads_per_sample = sub_df[samples].sum(axis=0)
            read_props = (taxon_reads_per_sample / total_reads_per_sample * 100).fillna(0).values
            res[taxon] = read_props
        res_df = pd.DataFrame(res, index=samples).transpose()

    else:
        # Case where specific categories from metadata are selected
        categories = metadata_df[selected_metadata].dropna().unique()

        # Filter metadata to keep relevant columns only
        metadata_df = metadata_df[['Sample', selected_metadata]]

        for taxon in all_taxa:
            sub_df = taxon_table_df[taxon_table_df[taxonomic_level] == taxon]

            reads_list = []
            for category in categories:
                # Get the samples belonging to the current category
                category_samples = metadata_df.loc[metadata_df[selected_metadata] == category, 'Sample'].tolist()

                # Calculate total and taxon-specific reads for the category
                total_reads = taxon_table_df[category_samples].sum().sum()
                taxon_reads = sub_df[category_samples].sum().sum()

                # Compute proportion
                read_props = (taxon_reads / total_reads * 100) if total_reads != 0 else 0
                reads_list.append(read_props)

            res[taxon] = reads_list
        res_df = pd.DataFrame(res, index=categories).transpose()

    # Use color map based on user settings
    colorscale_name = user_settings.get('colorsequence', 'Plotly')
    colors = getattr(px.colors.qualitative, colorscale_name, px.colors.qualitative.Plotly)

    # Initialize and populate the figure
    fig = go.Figure()

    for i, taxon in enumerate(res_df.index):
        fig.add_trace(go.Bar(
            x=res_df.columns,
            y=res_df.loc[taxon],
            name=taxon,
            hoverinfo='y+name',
            marker_color=colors[i % len(colors)]
        ))

    # Update layout
    fig.update_xaxes(dtick='linear')
    fig.update_yaxes(range=(0, 100))
    fig.update_layout(barmode='stack',
                      title=f'Proportion of Reads: {taxonomic_level}, {selected_metadata}',
                      template=user_settings['template'],
                      font_size=user_settings['font_size'],
                      yaxis_title='Proportion of Reads (%)',
                      yaxis_title_font=dict(size=user_settings['font_size']),
                      yaxis_tickfont=dict(size=user_settings['font_size']),
                      xaxis_title_font=dict(size=user_settings['font_size']),
                      xaxis_tickfont=dict(size=user_settings['font_size']))

    return fig


def read_props_heatmap(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, selected_metadata,
                       traits_df, selected_traits, user_settings, tool_settings):

    # Make copies of input dataframes to prevent altering the originals
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    # Extract relevant settings from the tool settings
    selected_metadata = tool_settings['selected_metadata']
    taxonomic_level = tool_settings['taxonomic_level']

    # Get unique taxa at the specified taxonomic level and exclude empty values
    all_taxa = taxon_table_df[taxonomic_level].dropna().unique()

    # Initialize result dictionary
    res = {}

    if selected_metadata == 'All samples':
        # Case where all samples are selected (simpler computation)
        total_reads_per_sample = taxon_table_df[samples].sum(axis=0)

        for taxon in all_taxa:
            # Filter taxon-specific data only once per taxon
            sub_df = taxon_table_df[taxon_table_df[taxonomic_level] == taxon]

            # Sum reads per taxon for each sample and calculate the proportion
            taxon_reads_per_sample = sub_df[samples].sum(axis=0)
            read_props = (taxon_reads_per_sample / total_reads_per_sample * 100).fillna(0).values
            res[taxon] = read_props
        res_df = pd.DataFrame(res, index=samples).transpose()

    else:
        # Case where specific categories from metadata are selected
        categories = [i for i in metadata_df[selected_metadata].dropna().unique() if i != '']

        # Filter metadata to keep relevant columns only
        metadata_df = metadata_df[['Sample', selected_metadata]]

        for taxon in all_taxa:
            sub_df = taxon_table_df[taxon_table_df[taxonomic_level] == taxon]

            reads_list = []
            for category in categories:
                # Get the samples belonging to the current category
                category_samples = metadata_df.loc[metadata_df[selected_metadata] == category, 'Sample'].tolist()

                # Calculate total and taxon-specific reads for the category
                total_reads = taxon_table_df[category_samples].sum().sum()
                taxon_reads = sub_df[category_samples].sum().sum()

                # Compute proportion
                read_props = (taxon_reads / total_reads * 100) if total_reads != 0 else 0
                reads_list.append(read_props)

            res[taxon] = reads_list
        res_df = pd.DataFrame(res, index=categories).transpose()

    # Use color map based on user settings
    colorscale_name = user_settings.get('colorscale', 'Viridis')

    # Initialize heatmap figure
    fig = go.Figure(data=go.Heatmap(
        z=res_df.values.tolist(),
        x=res_df.columns.tolist(),
        y=res_df.index.tolist(),
        hoverongaps=False,
        colorscale=colorscale_name))

    # Update layout
    fig.update_xaxes(dtick='linear')
    fig.update_yaxes(dtick='linear')
    fig.update_layout(
        title=f'Proportion of Reads: {taxonomic_level}, {selected_metadata}',
        template=user_settings['template'],
        font_size=user_settings['font_size'],
        yaxis_title='Proportion of Reads (%)',
        yaxis_title_font=dict(size=user_settings['font_size']),
        yaxis_tickfont=dict(size=user_settings['font_size']),
        xaxis_title_font=dict(size=user_settings['font_size']),
        xaxis_tickfont=dict(size=user_settings['font_size']))

    return fig













