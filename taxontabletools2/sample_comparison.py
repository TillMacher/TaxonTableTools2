import pandas as pd
import numpy as np
from pathlib import Path
from taxontabletools2.utilities import simple_taxontable
from taxontabletools2.utilities import filter_taxontable
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from io import BytesIO
import pandas as pd
import math
import numpy as np
from taxontabletools2.utilities import collect_traits
from taxontabletools2.utilities import strip_traits
from taxontabletools2.utilities import collect_metadata
from taxontabletools2.utilities import load_df
from taxontabletools2.utilities import collect_replicates
from taxontabletools2.utilities import export_taxon_table
from taxontabletools2.utilities import update_taxon_table
import requests, datetime, json, time
from stqdm import stqdm
import random, statistics
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def rarefaction_curve(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, selected_metadata, traits_df, selected_traits, users_settings, tool_settings):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    reps = tool_settings['reps']
    taxonomic_level = tool_settings['taxonomic_level']

    ## collect taxa per sample
    all_taxa = sorted(set(taxon_table_df[taxonomic_level].values.tolist()))
    taxa_dict = {}
    for sample in samples:
        hit_list = []
        for hit in taxon_table_df[[sample, taxonomic_level]].values.tolist():
            if hit[0] != 0 and hit[1] != '':
                hit_list.append(hit[1])
            taxa_dict[sample] = hit_list

    ## calculate read proportion of species
    reads_props_dict = {}
    total_reads = sum([sum(i) for i in taxon_table_df[samples].values.tolist()])
    for taxon in all_taxa:
        n_reads = taxon_table_df.loc[taxon_table_df[taxonomic_level] == taxon][samples].sum().sum() / total_reads * 100
        reads_props_dict[taxon] = n_reads

    ## perform rarefaction
    n_species_dict = {}
    n_reads_dict = {}
    stdev_dict = {}
    for j in range(1, len(samples)+1):
        n_species_list = []
        n_reads_list = []
        for i in range(1,reps):
            random_samples = random.sample(samples, k=j)
            species = set([item for sublist in [taxa_dict[sample] for sample in random_samples] for item in sublist])
            read_props = sum([reads_props_dict[i] for i in species])
            n_species = len(species)
            n_species_list.append(n_species)
            n_reads_list.append(read_props)

        stdev_dict[j] = statistics.stdev(n_species_list)
        n_species_dict[j] = np.mean(n_species_list)
        n_reads_dict[j] = np.mean(n_reads_list)

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    x_values = list(n_species_dict.keys())
    y_values = list(n_species_dict.values())
    y2_values = list(n_reads_dict.values())
    y_error = list(stdev_dict.values())

    ## add rarefaction
    fig.add_trace(go.Scatter(x=x_values, y=y_values, error_y=dict(type='data', array=y_error, visible=True, thickness=.5, width=.5), marker=dict(color='navy')))
    fig.add_trace(go.Scatter(x=x_values, y=y2_values, marker=dict(color='navy'), line_dash="dash"), secondary_y=True)

    ## calculate the point where less than one new species is found
    for t in [50, 75, 90, 100]:
        threshold = max(y_values) * t / 100
        cross_1 = [i for i in y_values if i >= threshold][0]
        x_1 = y_values.index(cross_1) + 1
        y_1 = cross_1
        if t == 100:
            tp = 'middle left'
        else:
            tp = 'middle right'
        fig.add_vline(x=x_1, line_width=1, line_dash="dash", line_color="grey")
        fig.add_trace(go.Scatter(x=[x_1], y=[y_1/2], mode='text', textfont=dict(size=7), textposition=tp, text=' {}% '.format(t, str(round(t,1))), line=dict(color='grey', width=1), name=str(t)))

    fig.update_layout(template='simple_white', title='', showlegend=False)
    fig.update_yaxes(rangemode="tozero", title=f'# {taxonomic_level.lower()}', secondary_y=False)
    fig.update_yaxes(range=(0,101), title='reads (%)', secondary_y=True)
    fig.update_xaxes(title='# samples', dtick=2, range=[1,len(samples)+1])

    out_df = pd.DataFrame()
    out_df['# samples'] = x_values
    out_df['# species'] = y_values
    out_df['reads (%)'] = y2_values

    return fig, out_df

def venn(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, selected_metadata, traits_df, selected_traits, users_settings, tool_settings):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    ## collect tool-specific settings
    selected_metadata = tool_settings['selected_metadata']
    taxonomic_level = tool_settings['taxonomic_level']

    ## collect the number of categories
    metadata_df = metadata_df[['Sample', selected_metadata]]
    categories = [i for i in metadata_df[selected_metadata].drop_duplicates().values.tolist() if i != '']
    n_categories = len(categories)

    if n_categories == 2:
        samples_a = metadata_df.loc[metadata_df[selected_metadata] == categories[0]]['Sample'].values.tolist()
        df_a = filter_taxontable(taxon_table_df, samples_a, taxonomic_level)
        species_a = set(df_a[taxonomic_level].values.tolist())
        samples_b = metadata_df.loc[metadata_df[selected_metadata] == categories[1]]['Sample'].values.tolist()
        df_b = filter_taxontable(taxon_table_df, samples_b, taxonomic_level)
        species_b = set(df_b[taxonomic_level].values.tolist())

        a_only = species_a - species_b
        print(a_only)
        n_a_only = len(a_only)
        shared = species_a & species_b
        print(shared)
        n_shared = len(shared)
        b_only = species_b - species_a
        print(b_only)
        n_b_only = len(b_only)

        # Create the Venn diagram
        plt.figure(figsize=(8, 8))
        venn = venn2(subsets=(n_a_only, n_b_only, n_shared), set_labels=(categories[0], categories[1]))

        # Adjust the font size of set labels
        venn.get_label_by_id('A').set_fontsize(12)  # Font size for "Category A"
        venn.get_label_by_id('B').set_fontsize(12)  # Font size for "Category B"
        # Adjust the font size of subset labels (the numbers inside the circles)
        for text in venn.subset_labels:
            if text:  # Check if the label exists (some subsets might be None)
                text.set_fontsize(users_settings['font_size'])

        return plt

    elif n_categories == 3:
        samples_a = metadata_df.loc[metadata_df[selected_metadata] == categories[0]]['Sample'].values.tolist()
        df_a = filter_taxontable(taxon_table_df, samples_a, taxonomic_level)
        species_a = set(df_a[taxonomic_level].values.tolist())

        samples_b = metadata_df.loc[metadata_df[selected_metadata] == categories[1]]['Sample'].values.tolist()
        df_b = filter_taxontable(taxon_table_df, samples_b, taxonomic_level)
        species_b = set(df_b[taxonomic_level].values.tolist())

        samples_c = metadata_df.loc[metadata_df[selected_metadata] == categories[2]]['Sample'].values.tolist()
        df_c = filter_taxontable(taxon_table_df, samples_c, taxonomic_level)
        species_c = set(df_c[taxonomic_level].values.tolist())

        # Calculate the different subsets
        a_only = species_a - (species_b | species_c)
        b_only = species_b - (species_a | species_c)
        c_only = species_c - (species_a | species_b)
        ab_shared = (species_a & species_b) - species_c
        ac_shared = (species_a & species_c) - species_b
        bc_shared = (species_b & species_c) - species_a
        abc_shared = species_a & species_b & species_c

        # Create the Venn diagram
        plt.figure(figsize=(8, 8))
        venn = venn3(subsets=(
        len(a_only), len(b_only), len(ab_shared), len(c_only), len(ac_shared), len(bc_shared), len(abc_shared)),
                     set_labels=(categories[0], categories[1], categories[2]))

        # Adjust the font size of set labels
        venn.get_label_by_id('A').set_fontsize(12)  # Font size for "Category A"
        venn.get_label_by_id('B').set_fontsize(12)  # Font size for "Category B"
        venn.get_label_by_id('C').set_fontsize(12)  # Font size for "Category C"

        # Adjust the font size of subset labels (the numbers inside the circles)
        for text in venn.subset_labels:
            if text:  # Check if the label exists (some subsets might be None)
                text.set_fontsize(users_settings['font_size'])

        return plt

    elif n_categories == 1:
        st.error(f'Venn diagrams only work with 2-3 categories. The metadata {selected_metadata} has one category.')
    else:
        st.error(
            f'Venn diagrams only work with 2-3 categories. The metadata {selected_metadata} has {n_categories} categories.')



