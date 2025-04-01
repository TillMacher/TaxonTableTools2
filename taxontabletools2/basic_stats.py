import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
import statistics

# TaXon_table_df = pd.read_excel('/Users/tillmacher/Desktop/TTT_projects/Projects/AA_test/TaXon_tables/tutorial_taxon_table.xlsx').fillna('nan')
# settings = {'path_to_projects': Path('/Users/tillmacher/Desktop/ttt_projects/Projects'), 'template': 'plotly', 'color1': 'blue', 'color2': 'black', 'colorsequence': 'Plotly', 'colorscale': 'jet', 'font_size': 16, 'clustering_unit': 'OTUs'}


@st.cache_data()
def basic_stats_reads(settings, TaXon_table_df):

    """ Calculate the number of reads per sample """
    samples = TaXon_table_df.columns.tolist()[10:]
    y_values = {i:sum([j for j in TaXon_table_df[i].values.tolist() if j != 0]) for i in samples}
    y_values = dict(sorted(y_values.items(), reverse=True, key=lambda item: item[1]))
    x_values = list(y_values.keys())

    max_reads, min_reads, avg_reads, stdev = max(y_values.values()), min(y_values.values()), round(statistics.mean(y_values.values()), 2), round(statistics.stdev(y_values.values()), 2)

    fig = go.Figure()
    fig = fig.add_trace(go.Bar(x=x_values,y=list(y_values.values())))
    fig.update_yaxes(title = 'Reads', title_font=dict(size=settings['font_size']), tickfont=dict(size=settings['font_size']))
    fig.update_xaxes(title='Samples', title_font=dict(size=settings['font_size']), showticklabels=False, tickfont=dict(size=settings['font_size']))
    fig.update_layout(template=settings['template'], font_size=settings['font_size'])

    return fig, min_reads, max_reads, avg_reads, stdev

@st.cache_data()
def basic_stats_OTUs(settings, TaXon_table_df):

    """ Calculate the number of OTUs per sample """

    samples = TaXon_table_df.columns.tolist()[10:]
    y_values = {i:len([j for j in TaXon_table_df[i].values.tolist() if j != 0]) for i in samples}
    y_values = dict(sorted(y_values.items(), reverse=True, key=lambda item: item[1]))
    x_values = list(y_values.keys())

    max_OTUs, min_OTUs, avg_OTUs = max(y_values.values()), min(y_values.values()), round(statistics.mean(y_values.values()), 2)

    fig = go.Figure()
    fig = fig.add_trace(go.Bar(x=x_values,y=list(y_values.values())))
    fig.update_layout(template=settings['template'], font_size=settings['font_size'])
    fig.update_yaxes(title = settings['clustering_unit'], title_font=dict(size=settings['font_size']), tickfont=dict(size=settings['font_size']))
    fig.update_xaxes(title='Samples', title_font=dict(size=settings['font_size']), showticklabels=False, tickfont=dict(size=settings['font_size']))

    return fig, min_OTUs, max_OTUs, avg_OTUs

@st.cache_data()
def taxonomic_richness(settings, TaXon_table_df):

    """ Calculate the taxonomic richness """

    taxonomic_levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    statistics_list, statistics_set, statistics_dict, highest_level_dict = [], [], {}, {}

    for taxon_to_evaluate in taxonomic_levels:
        taxa_list = [x for x in TaXon_table_df[taxon_to_evaluate].values.tolist() if str(x) != 'nan']
        statistics = taxon_to_evaluate, len(taxa_list)
        statistics_set.append(len(set(taxa_list)))
        statistics_list.append(list(statistics))
        statistics_dict[taxon_to_evaluate] = len(taxa_list)

    highest_level_dict["Phylum"] = statistics_dict["Phylum"] - statistics_dict["Class"]
    highest_level_dict["Class"] = statistics_dict["Class"] - statistics_dict["Order"]
    highest_level_dict["Order"] = statistics_dict["Order"] - statistics_dict["Family"]
    highest_level_dict["Family"] = statistics_dict["Family"] - statistics_dict["Genus"]
    highest_level_dict["Genus"] = statistics_dict["Genus"] - statistics_dict["Species"]
    highest_level_dict["Species"] = statistics_dict["Species"]

    taxon_levels = list(highest_level_dict.keys())
    number_of_taxa_per_level = statistics_set

    fig = go.Figure(data=[go.Bar(x=taxon_levels, y=number_of_taxa_per_level, name="Taxon", textposition="outside", cliponaxis=False, text=number_of_taxa_per_level)])
    fig.update_layout(template=settings['template'], font_size=settings['font_size'])
    fig.update_yaxes(title = '# taxa', title_font=dict(size=settings['font_size']), tickfont=dict(size=settings['font_size']))
    fig.update_xaxes(title_font=dict(size=settings['font_size']), tickfont=dict(size=settings['font_size']))

    return fig

@st.cache_data()
def taxonomic_resolution(settings, TaXon_table_df):
    """Calculate the taxonomic richness."""

    # All species
    all_species = TaXon_table_df['Species'].values.tolist()
    unique_species = sorted(set(all_species))
    species_count = {i:all_species.count(i) for i in unique_species if i != ''}
    species_count = dict(sorted(species_count.items(), key=lambda item: item[1], reverse=True))
    x_values = list(species_count.keys())[:15]
    display_x_values = []
    for value in x_values:
        try:
            display_x_values.append(f"<i>{value.split(' ')[0][0]}. {value.split(' ')[1]}<i>")
        except:
            display_x_values.append(value)

    y_values = list(species_count.values())[:15]

    # Title setup (not used in current code but retained for completeness)
    title = f"# {settings['clustering_unit']}"

    # Create bar chart
    fig = go.Figure(data=go.Bar(x=display_x_values, y=y_values))

    fig.update_layout(
        template=settings['template'],
        font_size=settings['font_size'],
        yaxis_title=settings['clustering_unit'],
    )

    fig.update_yaxes(title_font=dict(size=settings['font_size']), tickfont=dict(size=settings['font_size']))
    fig.update_xaxes(title_font=dict(size=settings['font_size']), tickfont=dict(size=settings['font_size']))

    return fig

#
