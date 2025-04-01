import pandas as pd
import math
import numpy as np
from pathlib import Path
import streamlit as st
from Bio.Seq import Seq
from Bio.Seq import CodonTable
from taxontabletools2.utilities import collect_traits
from taxontabletools2.utilities import strip_traits
from taxontabletools2.utilities import collect_metadata
from taxontabletools2.utilities import load_df
from taxontabletools2.utilities import collect_replicates
from taxontabletools2.utilities import export_taxon_table
from taxontabletools2.utilities import simple_taxontable
from taxontabletools2.utilities import update_taxon_table
from stqdm import stqdm

def presence_absence_table(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    ## set all element to 0/1
    taxon_table_df[samples] = taxon_table_df[samples].apply(lambda x: (x != 0).astype(int))

    ## export table
    suffix = 'PA'
    export_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, suffix)

def simplify_table(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df):
    # merge duplicate OTUs
    save = False
    filtered_df = simple_taxontable(taxon_table_xlsx, taxon_table_df, samples, metadata_df, save)
    ## export the dataframe
    new_TaXon_table_xlsx = str(Path(taxon_table_xlsx)).replace('.xlsx', '_simplified.xlsx')
    with pd.ExcelWriter(new_TaXon_table_xlsx) as writer:
        filtered_df.to_excel(writer, sheet_name='Taxon Table', index=False)
        metadata_df.to_excel(writer, sheet_name='Metadata Table', index=False)

def add_traits_from_file(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, new_traits_df, taxon_col):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    for row in new_traits_df.values.tolist():
        taxon = row[0]
        traits = row[1:]
        for trait, trait_name in zip(traits, new_traits_df.columns.tolist()[1:]):
            # Find the index of the rows in taxon_table_df that match the taxon
            indices = taxon_table_df[taxon_table_df[taxon_col] == taxon].index
            # Add traits to each matching row in taxon_table_df
            for index in indices:
                traits_df.loc[index, trait_name] = trait

    traits_df = traits_df.fillna('')
    update_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, '')

def sort_samples(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    # Collect the sorted samples
    sorted_samples = metadata_df['Sample'].values.tolist()
    taxon_table_df_sorted = taxon_table_df[taxon_table_df.columns.tolist()[:9] + sorted_samples]

    update_taxon_table(taxon_table_xlsx, taxon_table_df_sorted, traits_df, metadata_df, '')

def rename_samples(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, selected_metadata):
    # Create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    # Collect the new sample names from the selected_metadata column
    sorted_samples = metadata_df['Sample'].values.tolist()
    new_names = metadata_df[selected_metadata].values.tolist()

    # Rename columns in taxon_table_df using the old and new sample names
    for old, new in zip(sorted_samples, new_names):
        taxon_table_df.rename(columns={old: new}, inplace=True)

    # Update the metadata DataFrame with the old and new sample names
    metadata_df['Old names'] = sorted_samples
    metadata_df['Sample'] = new_names

    # Call the update function to apply changes
    update_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, '')

def merge_ESV_tables(TaXon_table_xlsx, df1, df2, df1_metadata, df2_metadata, suffix):
    df1_traits = collect_traits(df1)
    df1_trait_names = df1_traits.columns.tolist()[1:]
    df1 = strip_traits(df1)
    df1_samples = df1.columns.tolist()[9:]
    df2_traits = collect_traits(df2)
    df2_trait_names = df2_traits.columns.tolist()[1:]
    df2 = strip_traits(df2)
    df2_samples = df2.columns.tolist()[9:]

    ## calculate shared traits
    shared_traits = list(set(df1_trait_names) & set(df2_trait_names))

    if len(set(df1_samples) & set(df2_samples)) != 0:
        st.error('Error: Cannot merge tables with overlapping sample names!')

    else:
        ## Collect some stats
        ESVs_1 = set(df1['ID'].values.tolist())
        ESVs_2 = set(df2['ID'].values.tolist())
        a_only = len(ESVs_1 - ESVs_2)
        shared = len(ESVs_1 & ESVs_2)
        b_only = len(ESVs_2 - ESVs_1)
        total = a_only + shared + b_only
        all_IDs = sorted(set(df1['ID'].values.tolist() + df2['ID'].values.tolist()))

        st.success(f'Shared ESVs: {shared}')
        st.success(f'Original exclusive ESVs: {a_only}')
        st.success(f'Added exclusive ESVs: {b_only}')
        st.success(f'Total ESVs: {total}')

        merged_taxon_table_df = pd.DataFrame()
        for ID in all_IDs:
            new_row = pd.DataFrame()

            # Case 1: Present in both datasets: Merge them and select the taxonomy with the higher similarity
            if ID in ESVs_1 and ID in ESVs_2:
                # data from table 1
                row1 = df1.loc[df1['ID'] == ID]
                data1 = row1.iloc[:, :9]
                similarity1 = row1['Similarity'].values.tolist()[0]
                reads1 = row1[df1_samples]
                traits1 = df1_traits.loc[df1_traits['ID'] == ID][shared_traits]
                traits1['ESV merge'] = ['Shared']

                # data from table 2
                row2 = df2.loc[df2['ID'] == ID]
                data2 = row2.iloc[:, :9]
                similarity2 = row2['Similarity'].values.tolist()[0]
                reads2 = row2[df2_samples]
                traits2 = df2_traits.loc[df2_traits['ID'] == ID][shared_traits]

                if data1.values.tolist() == data2.values.tolist():
                    # Find the index of the column "Seq"
                    col_index = data1.columns.get_loc("Seq")
                    # Insert the new DataFrame columns before "Seq"
                    data1_traits = pd.concat([data1.iloc[:, :col_index], traits1, data1.iloc[:, col_index:]], axis=1)
                    new_row = pd.concat([data1_traits.reset_index(drop=True),
                                         reads1.reset_index(drop=True),
                                         reads2.reset_index(drop=True)], axis=1)

                elif similarity1 > similarity2:
                    # Find the index of the column "Seq"
                    col_index = data1.columns.get_loc("Seq")
                    # Insert the new DataFrame columns before "Seq"
                    data1_traits = pd.concat([data1.iloc[:, :col_index], traits1, data1.iloc[:, col_index:]], axis=1)
                    new_row = pd.concat([data1_traits.reset_index(drop=True),
                                         reads1.reset_index(drop=True),
                                         reads2.reset_index(drop=True)], axis=1)
                else:
                    # Find the index of the column "Seq"
                    col_index = data2.columns.get_loc("Seq")
                    # Insert the new DataFrame columns before "Seq"
                    data2_traits = pd.concat([data2.iloc[:, :col_index], traits2, data2.iloc[:, col_index:]], axis=1)
                    new_row = pd.concat([data2_traits.reset_index(drop=True),
                                         reads1.reset_index(drop=True),
                                         reads2.reset_index(drop=True)], axis=1)

            # Case 2: Present in dataset 1: Fill up the other
            elif ID in ESVs_1:
                # data from table 1
                row1 = df1.loc[df1['ID'] == ID]
                data1 = row1.iloc[:, :9]
                reads1 = row1[df1_samples]
                traits1 = df1_traits.loc[df1_traits['ID'] == ID][shared_traits]
                traits1['ESV merge'] = ['Original']

                # data from table 2
                reads2 = pd.DataFrame([[0] * len(df2_samples)], columns=df2_samples)

                # Find the index of the column "Seq"
                col_index = data1.columns.get_loc("Seq")
                # Insert the new DataFrame columns before "Seq"
                data1_traits = pd.concat([data1.iloc[:, :col_index], traits1, data1.iloc[:, col_index:]], axis=1)
                new_row = pd.concat([data1_traits.reset_index(drop=True),
                                     reads1.reset_index(drop=True),
                                     reads2.reset_index(drop=True)], axis=1)

            # Case 2: Present in dataset 1: Fill up the other
            elif ID in ESVs_2:
                # data from table 1
                row2 = df2.loc[df2['ID'] == ID]
                data2 = row2.iloc[:, :9]
                reads2 = row2[df2_samples]
                traits2 = df2_traits.loc[df2_traits['ID'] == ID][shared_traits]
                traits2['ESV merge'] = ['Added']

                # data from table 2
                reads1 = pd.DataFrame([[0] * len(df1_samples)], columns=df1_samples)

                # Find the index of the column "Seq"
                col_index = data2.columns.get_loc("Seq")
                # Insert the new DataFrame columns before "Seq"
                data2_traits = pd.concat([data2.iloc[:, :col_index], traits2, data2.iloc[:, col_index:]], axis=1)
                new_row = pd.concat([data2_traits.reset_index(drop=True),
                                     reads1.reset_index(drop=True),
                                     reads2.reset_index(drop=True)], axis=1)

            else:
                print('Whoops that should not happen!')

            ## add to newly merged df
            merged_taxon_table_df = pd.concat([merged_taxon_table_df, new_row], ignore_index=True)

    ## the merged_df is the new Taxon Table
    ## now create the Metadata Table
    ## calculate shared metadata
    shared_metadata = list(set(df1_metadata.columns.tolist()) & set(df2_metadata.columns.tolist()))
    export_metadata_table_df = pd.concat([df1_metadata[shared_metadata], df2_metadata[shared_metadata]], ignore_index=True)
    export_metadata_table_df['ESV merge'] = ['Original'] * len(df1_metadata) + ['Added'] * len(df2_metadata)

    ## export df
    export_traits_df = collect_traits(merged_taxon_table_df)
    export_taxon_table_df = strip_traits(merged_taxon_table_df)
    export_taxon_table(TaXon_table_xlsx, export_taxon_table_df, export_traits_df, export_metadata_table_df, suffix)

    st.success(f'Tables were successfully merged!')

def export_fasta(taxon_table_xlsx, taxon_table_df, user_settings):
    fasta_file = user_settings['path_to_outdirs'].joinpath(user_settings['project_name'], 'Fasta', taxon_table_xlsx.stem + '.fasta')

def open_reading_frame(taxon_table_xlsx, taxon_table_df, metadata_df, traits_df, translation_table, orf_filter):

    def multiple_of_three(dna_seq):
        for i in range(0, 3):
            c = (len(dna_seq) + i) % 3
            if c == 0:
                return dna_seq + Seq('N' * i)

    sequences = taxon_table_df['Seq'].values.tolist()
    ids = taxon_table_df['ID'].values.tolist()

    frames = []
    for id, sequence in stqdm(zip(ids, sequences)):
        sequence = sequence.upper()
        if all(c in 'ACGTN' for c in sequence):
            dna_seq = Seq(sequence)
            res = [id]
            stops = []
            for i in range(3):
                seq = multiple_of_three(dna_seq[i:])
                prot = seq.translate(table=translation_table, to_stop=False)
                n_stop = str(prot).count('*')
                stops.append(n_stop)
                res.extend([n_stop, str(prot)])

            # find best reading frame
            pos = stops.index(min(stops))
            if min(stops) == 0:
                res.extend([pos+1])
            else:
                res.extend([''])
            frames.append(res)
        else:
            frames.append([id] + ['']*7)

    # Create dataframe
    frames_df = pd.DataFrame(frames, columns=['ID', 'stop1', 'prot1', 'stop2', 'prot2', 'stop3', 'prot3', 'clean frame'])
    # Remove old columns
    if set(['stop1', 'prot1', 'stop2', 'prot2', 'stop3', 'prot3', 'clean frame']).issubset(traits_df.columns):
        traits_df.drop(['stop1', 'prot1', 'stop2', 'prot2', 'stop3', 'prot3', 'clean frame'], axis=1, inplace=True)
    # Merge with traits
    traits_df = pd.merge(frames_df, traits_df, on='ID', how='inner')

    if orf_filter == 'Yes':
        traits_df_filtered = traits_df.loc[traits_df['clean frame'] != '']
        taxon_table_df_filtered = taxon_table_df.loc[taxon_table_df['ID'].isin(traits_df_filtered['ID'])]
        export_taxon_table(taxon_table_xlsx, taxon_table_df_filtered, traits_df_filtered, metadata_df, 'orf')
        st.success(f'Filtered tables for sequences with open reading frames! Removed {len(taxon_table_df) - len(taxon_table_df_filtered)} sequences.')
    else:
        update_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, '')
        st.success('Added reading frames to taxon table!')













