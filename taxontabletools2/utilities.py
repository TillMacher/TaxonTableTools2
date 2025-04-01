import streamlit as st
import pandas as pd
import numpy as np
import scipy, subprocess, webbrowser
from pathlib import Path
import glob, sys, time, statistics, os.path
import seqconverter


########################################################################################################################
def load_df(table):
    df = pd.read_excel(table).fillna('')
    return df

def reduce_taxontable(df, taxonomic_level, samples, meta_data_to_test):

    if meta_data_to_test == '':
        ## extract the relevant data
        df = df[[taxonomic_level] + samples]
        ## define an aggregation function to combine multiple hit of one taxonimic level
        aggregation_functions = {}
        ## define samples functions
        for sample in samples:
            ## 'sum' will calculate the sum of p/a data
            aggregation_functions[sample] = 'sum'
        ## define taxon level function
        aggregation_functions[taxonomic_level] = 'first'
        ## create condensed dataframe
        df = df.groupby(df[taxonomic_level]).aggregate(aggregation_functions)
        if '' in df.index:
            df = df.drop('')

        return df

    else:
        ## collect metadata for each taxon
        taxon_metadata_dict = {}
        for row in df[[taxonomic_level, meta_data_to_test]].values.tolist():
            taxon = row[0]
            metadata = row[1]
            taxon_metadata_dict[taxon] = metadata

        ## extract the relevant data
        df = df[[taxonomic_level] + samples]
        ## define an aggregation function to combine multiple hit of one taxonimic level
        aggregation_functions = {}
        ## define samples functions
        for sample in samples:
            ## 'sum' will calculate the sum of p/a data
            aggregation_functions[sample] = 'sum'
        ## define taxon level function
        aggregation_functions[taxonomic_level] = 'first'
        ## create condensed dataframe
        df = df.groupby(df[taxonomic_level]).aggregate(aggregation_functions)
        if '' in df.index:
            df = df.drop('')

        ## add metadata back to taxa
        metadata_col = []
        for taxon in df[taxonomic_level].values.tolist():
            metadata_col.append(taxon_metadata_dict[taxon])
        df[meta_data_to_test] = metadata_col

        return df

def strip_traits(df):
    " Strip the additional metadata from the dataframe to produce a clean TaxonTable "

    ## position of seq == acts as separator
    seq_loc = df.columns.get_loc("Seq")

    ## samples are always assending the seq column
    samples = df.columns.tolist()[seq_loc+1:]

    ## standard columns
    standard_columns = ['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'Seq']

    ## extract legacy TTT frame
    df = df[standard_columns + samples]

    return df

def collect_traits(df):
    " Collect the additional metadata from the dataframe to be able to place it back later "

    ## position of status and seq == acts as separator
    seq_loc = df.columns.get_loc("Seq")
    status = df.columns.get_loc("Similarity")

    ## metadata is always located between status and seq
    metadata = df.columns.tolist()[status+1:seq_loc]

    ## collect IDs and metadata
    df = df[['ID'] + metadata]

    return df

def add_traits(df, traits_df):
    " Collect the additional metadata from the dataframe to add them back to a df "

    for metadata in  traits_df.columns.tolist()[1:]:
        status_loc = df.columns.get_loc("Seq")
        ## loope through IDs and find all that still remain in the altered Taxon table (e.g. after filtering less OTUs might be present)
        sub_df = pd.DataFrame([i for i in traits_df[['ID', metadata]].values.tolist() if i[0] in df['ID'].values.tolist()], columns = ['ID', metadata])
        ## insert sub_df containing the respective metadata after the Status column
        df.insert(status_loc, str(metadata), sub_df[metadata].values.tolist(), True)

    return df

def collect_metadata(taxon_table_xlsx):
    try:
        metadata_df = pd.read_excel(taxon_table_xlsx, sheet_name='Metadata table')
    except ValueError:
        metadata_df = pd.DataFrame(columns=['Sample', 'Metadata'])
    return metadata_df

def export_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, suffix):
    ## collect remaining samples and filter metadata_df
    seq_loc = taxon_table_df.columns.tolist().index('Seq')
    samples = taxon_table_df.columns.tolist()[seq_loc+1::]
    filtered_metadata_df = metadata_df.loc[metadata_df['Sample'].isin(samples)]

    ## merge taxon table and trait table
    merged_taxon_table_df = taxon_table_df.merge(traits_df, on='ID', how='inner')
    cols = list(merged_taxon_table_df.columns)
    seq_index = cols.index('Seq')
    traits_cols = traits_df.columns.tolist()
    traits_cols.remove('ID')  # remove 'ID' as it is already in `taxon_table_df`
    [cols.remove(trait) for trait in traits_cols]
    new_cols = cols[:seq_index] + traits_cols + cols[seq_index:]
    merged_taxon_table_df = merged_taxon_table_df[new_cols]

    ## Calculate the number of empty OTUs
    empty_otus_count = (~(merged_taxon_table_df[samples] != 0).any(axis=1)).sum()

    ## Remove empty OTUs from dataframe
    mask = (merged_taxon_table_df[samples] != 0).any(axis=1)
    merged_taxon_table_df = merged_taxon_table_df[mask]

    if empty_otus_count != 0:
        st.warning(f'Warning: Removed {empty_otus_count} IDs where all samples were empty (zero reads).')

    ## create a new file name
    new_TaXon_table_xlsx = Path(str(Path(taxon_table_xlsx)).replace('.xlsx', f'_{suffix}.xlsx'))

    ## sort taxon table
    merged_taxon_table_df = merged_taxon_table_df.sort_values(by=['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity'])

    ## export the dataframe
    with pd.ExcelWriter(new_TaXon_table_xlsx) as writer:
        merged_taxon_table_df.to_excel(writer, sheet_name='Taxon Table', index=False)
        filtered_metadata_df.to_excel(writer, sheet_name='Metadata Table', index=False)

def update_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, suffix):
    ## collect remaining samples and filter metadata_df
    seq_loc = taxon_table_df.columns.tolist().index('Seq')
    samples = taxon_table_df.columns.tolist()[seq_loc+1::]
    filtered_metadata_df = metadata_df.loc[metadata_df['Sample'].isin(samples)]

    ## merge taxon table and trait table
    merged_taxon_table_df = taxon_table_df.merge(traits_df, on='ID', how='inner')
    cols = list(merged_taxon_table_df.columns)
    seq_index = cols.index('Seq')
    traits_cols = traits_df.columns.tolist()
    traits_cols.remove('ID')  # remove 'ID' as it is already in `taxon_table_df`
    [cols.remove(trait) for trait in traits_cols]
    new_cols = cols[:seq_index] + traits_cols + cols[seq_index:]
    merged_taxon_table_df = merged_taxon_table_df[new_cols]

    ## remove empty OTUs from dataframe
    mask = (merged_taxon_table_df[samples] != 0).any(axis=1)
    merged_taxon_table_df = merged_taxon_table_df[mask]

    ## export the dataframe
    with pd.ExcelWriter(taxon_table_xlsx) as writer:
        merged_taxon_table_df.to_excel(writer, sheet_name='Taxon Table', index=False)
        filtered_metadata_df.to_excel(writer, sheet_name='Metadata Table', index=False)

def simple_taxontable(taxon_table_xlsx, taxon_table_df, samples, metadata_df, save):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()

    ## add unique taxa
    taxon_table_df['Taxon'] = [' '.join(i) for i in taxon_table_df[['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].values.tolist()]
    unique_taxa = taxon_table_df['Taxon'].drop_duplicates().values.tolist()
    new_df_list = []
    n_OTUs_list = []
    for i, taxon in enumerate(unique_taxa):
        sub_df = taxon_table_df.loc[taxon_table_df['Taxon'] == taxon].copy()
        sub_df.pop('Taxon')
        sub_df.loc[:, 'sub_df_sort'] = sub_df[samples].sum(axis=1)
        sub_df = sub_df.sort_values('sub_df_sort', ascending=False)
        sub_df.pop('sub_df_sort')
        n_OTUs_list.append(len(sub_df))
        samples = sub_df.columns.tolist()[sub_df.columns.tolist().index('Seq')+1:]
        reads = sub_df[samples].sum(axis=0)
        OTU_name = f'spOTU_{i}'
        new_row = [OTU_name] + sub_df.iloc[0].values.tolist()[1:sub_df.columns.tolist().index('Seq')+1] + reads.values.tolist()
        new_df_list.append(new_row)
    taxon_table_df.pop('Taxon')
    new_df = pd.DataFrame(new_df_list, columns=taxon_table_df.columns.tolist())
    new_df['Similarity'] = 100
    new_df.insert(9, "Merged", n_OTUs_list)
    new_df['Seq'] = 'Please refer to the non-simplified table!'

    ## create a new file name
    new_TaXon_table_xlsx = str(Path(taxon_table_xlsx)).replace('.xlsx', '_simplified.xlsx')

    if save == True:
        ## export the dataframe
        with pd.ExcelWriter(new_TaXon_table_xlsx) as writer:
            new_df.to_excel(writer, sheet_name='Taxon Table', index=False)
            metadata_df.to_excel(writer, sheet_name='Metadata Table', index=False)
    else:
        return new_df

def filter_taxontable(taxon_table_df, sample_subset, taxonomic_level):
    sub_df = taxon_table_df[[taxonomic_level] + sample_subset]
    sub_df = sub_df.loc[sub_df[taxonomic_level] != '']
    merged_df = sub_df.groupby(taxonomic_level, as_index=False).sum()
    filtered_df = merged_df[(merged_df.iloc[:, 1:].sum(axis=1) != 0)]
    return  filtered_df

def filter_taxontable_2(taxon_table_df, sample_subset, taxonomic_level_1, taxonomic_level_2):
    sub_df = taxon_table_df[[taxonomic_level_1, taxonomic_level_2] + sample_subset]
    sub_df = sub_df[(sub_df[taxonomic_level_1] != '') | (sub_df[taxonomic_level_2] != '')]
    sub_df['Combined_Taxon'] = sub_df[taxonomic_level_1] + ' ' + sub_df[taxonomic_level_2]
    merged_df = sub_df.groupby(['Combined_Taxon', taxonomic_level_1, taxonomic_level_2], as_index=False).sum()
    filtered_df = merged_df[(merged_df.iloc[:, 3:].sum(axis=1) != 0)]
    filtered_df.drop(columns=['Combined_Taxon'], inplace=True)
    column_order = [taxonomic_level_1, taxonomic_level_2] + sample_subset
    filtered_df = filtered_df[column_order]
    return filtered_df

########################################################################################################################

def taxon_richness_df(taxon_table_df, samples, taxonomic_level_1):
    richness_dict = {}  # Number of taxa per sample
    for sample in samples:
        sub_df = taxon_table_df[[sample, taxonomic_level_1]]
        sub_df = sub_df.loc[sub_df[sample] != 0]
        richness_dict[sample] = [len(sub_df[taxonomic_level_1].drop_duplicates())]

    richness_df = pd.DataFrame(richness_dict, index=[taxonomic_level_1])
    return richness_df

def taxon_reads_df(taxon_table_df, samples, taxonomic_level_1):
    reads_dict = {}  # Number of reads per taxon per sample
    all_taxa = sorted([i for i in taxon_table_df[taxonomic_level_1].drop_duplicates().values.tolist() if i != ''])

    for sample in samples:
        sub_df = taxon_table_df[[sample, taxonomic_level_1]]
        sub_df = sub_df.loc[sub_df[sample] != 0]

        reads_list = []
        for taxon in all_taxa:
            taxon_df = sub_df.loc[sub_df[taxonomic_level_1] == taxon]
            reads = taxon_df[sample].sum()
            reads_list.append(reads)

        reads_dict[sample] = reads_list

    reads_df = pd.DataFrame(reads_dict, index=all_taxa)
    return reads_df

def taxon_read_proportions_df(taxon_table_df, samples, taxonomic_level_1):
    read_proportions_dict = {}  # Read proportion per taxon per sample
    all_taxa = sorted([i for i in taxon_table_df[taxonomic_level_1].drop_duplicates().values.tolist() if i != ''])

    for sample in samples:
        sub_df = taxon_table_df[[sample, taxonomic_level_1]]
        sub_df = sub_df.loc[sub_df[sample] != 0]
        sample_reads = sub_df[sample].sum()

        read_proportions_list = []
        for taxon in all_taxa:
            taxon_df = sub_df.loc[sub_df[taxonomic_level_1] == taxon]
            reads = taxon_df[sample].sum()
            read_proportions_list.append(reads / sample_reads * 100 if sample_reads != 0 else 0)

        read_proportions_dict[sample] = read_proportions_list

    read_proportions_df = pd.DataFrame(read_proportions_dict, index=all_taxa)
    return read_proportions_df

def taxon_occurrence_df(taxon_table_df, samples, taxonomic_level_1):
    occurrence_dict = {}  # Occurrence (1/0) per taxon per sample
    all_taxa = sorted([i for i in taxon_table_df[taxonomic_level_1].drop_duplicates().values.tolist() if i != ''])

    for sample in samples:
        sub_df = taxon_table_df[[sample, taxonomic_level_1]]
        sub_df = sub_df.loc[sub_df[sample] != 0]

        occurrence_list = []
        for taxon in all_taxa:
            taxon_df = sub_df.loc[sub_df[taxonomic_level_1] == taxon]
            reads = taxon_df[sample].sum()
            occurrence_list.append(1 if reads != 0 else 0)

        occurrence_dict[sample] = occurrence_list

    occurrence_df = pd.DataFrame(occurrence_dict, index=all_taxa)
    return occurrence_df

def taxon_richness_per_higher_taxon_df(taxon_table_df, samples, taxonomic_level_1, taxonomic_level_2):
    richness_per_taxon_dict = {}  # Richness per higher taxon (e.g., family)
    all_taxa = sorted([i for i in taxon_table_df[taxonomic_level_2].drop_duplicates().values.tolist() if i != ''])

    for sample in samples:
        sub_df = taxon_table_df[[sample, taxonomic_level_1, taxonomic_level_2]]
        sub_df = sub_df.loc[sub_df[sample] != 0]

        richness_per_taxon_list = []
        for taxon in all_taxa:
            taxon_df = sub_df.loc[sub_df[taxonomic_level_2] == taxon]
            richness_per_taxon_list.append(len(taxon_df[taxonomic_level_1].drop_duplicates()))

        richness_per_taxon_dict[sample] = richness_per_taxon_list

    richness_per_taxon_df = pd.DataFrame(richness_per_taxon_dict, index=all_taxa)
    return richness_per_taxon_df if taxonomic_level_2 else None

########################################################################################################################
def collect_replicates(samples):

    unique_names = sorted(list(set(['_'.join(i.split('_')[:-1]) for i in samples])))
    suffixes = sorted(list(set([i.split('_')[-1] for i in samples])))

    return [unique_names, suffixes]

def check_replicates(replicates_dict, samples):
    unique_samples, replicate_suffixes = replicates_dict
    expected_samples = [sample+'_'+suffix for sample in unique_samples for suffix in replicate_suffixes]
    missing_samples = set(expected_samples) - set(samples)
    test = sorted(samples) == sorted(expected_samples)
    return test, missing_samples

########################################################################################################################

# taxon_table_xlsx = Path('/Users/tillmacher/Desktop/TTT_projects/AA_test/TaXon_tables/Ruwer_eDNA_fish_2024_d4_taxon_table_RuMo_merged_NCsub_fish_norm.xlsx')
# folder = "Venn_diagrams"
# file_name = "PCoA"
# fig = ""
# settings = {'project_name': 'AA_test', 'path_to_outdirs': Path('/Users/tillmacher/Desktop/TTT_projects'), 'plot_height': 1000, 'plot_width': 1000, 'show_legend': True, 'template': 'simple_white', 'font_size': 20, 'clustering_unit': 'ESVs', 'scatter_size': 15, 'color_1': 'navy', 'color_2': 'teal', 'colorsequence': 'Plotly', 'colorscale': 'blues'}
# lib = "matplot"

def export_plot(taxon_table_xlsx, folder, file_name, suffix, fig, settings, lib):
    xlsx_name = Path(taxon_table_xlsx).stem
    path_to_outdirs = settings['path_to_outdirs']

    if lib == 'matplot':
        file_pdf = Path(str(path_to_outdirs)).joinpath(settings['project_name'], folder, f'{xlsx_name}_{suffix}_{file_name}.pdf')
        fig.savefig(file_pdf, dpi=300, bbox_inches="tight")
        st.success(f'Saved plot as pdf!')

    if lib == 'plotly':
        file_pdf = Path(str(path_to_outdirs)).joinpath(settings['project_name'], folder, f'{xlsx_name}_{suffix}_{file_name}.pdf')
        file_html = Path(str(path_to_outdirs)).joinpath(settings['project_name'], folder, f'{xlsx_name}_{suffix}_{file_name}.html')
        fig.update_layout(height=settings['plot_height'], width=settings['plot_width'])
        fig.write_image(file_pdf)
        fig.write_html(file_html)
        st.success(f'Saved plot as pdf and html!')

########################################################################################################################
def spearman(a, b):
    # spearman's rho
    results = scipy.stats.spearmanr(a, b)
    spearman_p = results[1]
    if spearman_p <= 0.05:
        spearman_rho = str(round(results[0], 3)) + "*"
    else:
        spearman_rho = str(round(results[0], 3))
    return spearman_rho

########################################################################################################################
def glob_available_taxon_tables(path_to_outdirs):
    TaXon_tables_dict = {}
    tables = sorted(glob.glob(str(Path(path_to_outdirs).joinpath('TaXon_tables', '*.xlsx'))), key=os.path.getmtime)[::-1]
    for table in tables:
        TaXon_tables_dict[Path(table).stem] = Path(table)
    return TaXon_tables_dict

########################################################################################################################
def collect_sample_stats(TaXon_table_df, TaXon_table_name, traits_df, samples, metadata_df):
    ## information about the table
    n_samples = len(samples)
    n_OTUs = len(TaXon_table_df['ID'])
    n_phyla = len(set([i for i in TaXon_table_df['Phylum'].values.tolist() if i != '']))
    n_classes = len(set([i for i in TaXon_table_df['Class'].values.tolist() if i != '']))
    n_orders = len(set([i for i in TaXon_table_df['Order'].values.tolist() if i != '']))
    n_families = len(set([i for i in TaXon_table_df['Family'].values.tolist() if i != '']))
    n_genera = len(set([i for i in TaXon_table_df['Genus'].values.tolist() if i != '']))
    n_species = len(set([i for i in TaXon_table_df['Species'].values.tolist() if i != '']))
    n_traits = len(traits_df.columns[1:])
    n_reads = sum([sum(i[10:]) for i in TaXon_table_df.values.tolist()])
    cols = ['Samples', 'ID', 'Phyla', 'Classes', 'Orders', 'Families', 'Genera', 'Species', 'Traits',
            'Reads']
    row = [[n_samples, n_OTUs, n_phyla, n_classes, n_orders, n_families, n_genera, n_species, n_traits, n_reads]]
    stats_df = pd.DataFrame(row, columns=cols, index=['#'])

    # collect information about each sample in the table
    samples_info_list = []
    for sample in samples:
        tmp = pd.DataFrame([i for i in TaXon_table_df[
            ['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', sample]].values.tolist() if i[-1] != 0],
                           columns=['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', sample])
        n_OTUs = len(tmp['ID'])
        n_phyla = len(set([i for i in tmp['Phylum'].values.tolist() if i != '']))
        n_classes = len(set([i for i in tmp['Class'].values.tolist() if i != '']))
        n_orders = len(set([i for i in tmp['Order'].values.tolist() if i != '']))
        n_families = len(set([i for i in tmp['Family'].values.tolist() if i != '']))
        n_genera = len(set([i for i in tmp['Genus'].values.tolist() if i != '']))
        n_species = len(set([i for i in tmp['Species'].values.tolist() if i != '']))
        n_reads = sum(tmp[sample])
        samples_info_list.append(
            [sample, n_OTUs, n_phyla, n_classes, n_orders, n_families, n_genera, n_species, n_reads])

    sample_stats_df = pd.DataFrame(samples_info_list,
                                   columns=['Samples', 'ID', 'Phyla', 'Classes', 'Orders',
                                            'Families', 'Genera', 'Species', 'Reads'])

    # Calculate the mean for each numeric column
    averages = sample_stats_df.select_dtypes(include=[np.number]).mean()
    # Append the averages to the DataFrame using pandas.concat
    sample_stats_df = pd.concat([pd.DataFrame(averages).T, sample_stats_df], ignore_index=True)
    # Set the 'Samples' column of the first row to 'Average'
    sample_stats_df.loc[0, 'Samples'] = 'Average'
    # Replace the index with the 'Samples' column
    sample_stats_df.set_index('Samples', inplace=True)

    # Function to format numbers
    def format_numbers(row):
        if row.name == 'Average':
            return row.round(2).astype(str)
        else:
            return row.astype(int).astype(str)

    # Apply the function to each row
    sample_stats_df = sample_stats_df.apply(format_numbers, axis=1)

    return sample_stats_df

def taxon_table_fasta(taxon_table_xlsx, taxon_table_df):
    taxon_table_xlsx = '/Users/tillmacher/Desktop/TTT_projects/Projects/AA_test/TaXon_tables/All_ESVs_ID_merged_ID_merge_ID_merge_ID_merge_norm.xlsx'
    taxon_table_df = pd.read_excel('/Users/tillmacher/Desktop/TTT_projects/Projects/AA_test/TaXon_tables/All_ESVs_ID_merged_ID_merge_ID_merge_ID_merge_norm.xlsx').fillna('')
    taxon_table_fasta = taxon_table_xlsx.replace('.xlsx', '.fasta')
    f = open(taxon_table_fasta, 'w')
    for row in taxon_table_df[['ID', 'Seq']].values.tolist():
        f.write(f'>{row[0]}\n')
        f.write(f'{row[1]}\n')
    f.close()

########################################################################################################################








