import pandas as pd
import math
import numpy as np
from taxontabletools2.utilities import collect_traits
from taxontabletools2.utilities import strip_traits
from taxontabletools2.utilities import collect_metadata
from taxontabletools2.utilities import load_df
from taxontabletools2.utilities import collect_replicates
from taxontabletools2.utilities import export_taxon_table

def sum_if_less_than_X_zeros(row, cutoff):
    if (row != 0).sum() < cutoff:
        return 0
    else:
        return row.sum()

def replicate_merging(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, cutoff, replicates_dict):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    ## collect clean sample names
    unique_samples, suffixes = replicates_dict
    clean_samples = sorted(set([sample.replace('_' + suffix, '') for sample in samples for suffix in suffixes if suffix in sample]))

    ## merge replicates
    for sample in clean_samples:
        replicates = [f'{sample}_{i}' for i in suffixes]
        if all(value in samples for value in replicates):
            sub_df = taxon_table_df[replicates]
            taxon_table_df[sample] = sub_df.apply(sum_if_less_than_X_zeros, axis=1, args=(cutoff,))
            taxon_table_df.drop(columns=replicates, inplace=True)

    suffix = 'merged'
    export_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, suffix)

def negative_control_subtraction(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, mode, NC_list):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()


    ## collect samples
    samples_keep = [i for i in samples if i not in NC_list]
    # filter for samples
    sample_df = taxon_table_df[samples_keep]

    ## calculate NC reads
    if mode == 'Sum of all NCs':
        NC_df = taxon_table_df[NC_list].sum(axis=1)
    elif mode == 'Maximum in all NCs':
        NC_df = taxon_table_df[NC_list].max(axis=1)
    elif mode == 'Average in all NCs':
        NC_df = taxon_table_df[NC_list].mean(axis=1)
    result_df = sample_df.subtract(NC_df, axis=0)
    result_df = result_df.clip(lower=0)
    result_df['ID'] = taxon_table_df['ID']

    ## merge tables again
    taxon_table_df.drop(samples, axis=1, inplace=True)
    merged_taxon_table_df = taxon_table_df.merge(result_df, on='ID', how='inner')

    ## export table
    suffix = 'NCsub'
    export_taxon_table(taxon_table_xlsx, merged_taxon_table_df, traits_df, metadata_df, suffix)

def read_based_filter(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, mode, read_cutoff):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()


    ## filter reads
    if mode == 'Absolute':
        taxon_table_df[samples] = taxon_table_df[samples].where(taxon_table_df[samples] >= read_cutoff, 0)
    elif mode == 'Relative':
        for sample in samples:
            sample_cutoff = math.ceil(taxon_table_df[sample].sum() * read_cutoff)
            taxon_table_df[sample] = taxon_table_df[sample].where(taxon_table_df[sample] >= sample_cutoff, 0)

    ## export table
    suffix = 'filtered'
    export_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, suffix)

def read_based_normalisation(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, sub_sample_size):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    ## count the lowest number of reads to use as threshold
    reads_list = [sum(taxon_table_df[i].values.tolist()) for i in samples]

    ## convert sub_sample_size to int
    sub_sample_size = int(sub_sample_size)

    ## collect all OTUs
    OTU_list = taxon_table_df["ID"].values.tolist()

    df_out = pd.DataFrame()
    required_columns = taxon_table_df.columns.tolist()[0:10]
    df_out[required_columns] = taxon_table_df[required_columns]

    warnings = []

    for sample in samples:
        ## filter sample from data
        read_df = taxon_table_df[[sample, "ID"]]
        ## drop empty OTUs
        read_df = read_df[read_df[sample] != 0]
        ## check if sample can be normalized, otherwise just keep all reads and OTUs
        if sub_sample_size <= sum(read_df[sample].values.tolist()):
            ## create read list to draw the subsamples from
            read_list = pd.Series(np.repeat(read_df['ID'].to_list(), read_df[sample].to_list()))
            ## create a random subsample
            sub_sample = read_list.sample(n = sub_sample_size)
            ## count the number of reads per OTU
            sub_sample_reads = dict(pd.Series(sub_sample).value_counts())
            ## create a sorted OTU list
            OTU_sample_list = []
            for OTU in OTU_list:
                if OTU in sub_sample_reads.keys():
                    OTU_sample_list.append(sub_sample_reads[OTU])
                else:
                    OTU_sample_list.append(0)
        else:
            OTU_sample_list = taxon_table_df[[sample]]
            warnings.append(sample)

        ## add OTUs to dataframe
        df_out[sample] = OTU_sample_list

    ## export table
    suffix = 'norm'
    export_taxon_table(taxon_table_xlsx, df_out, traits_df, metadata_df, suffix)

def taxonomic_filtering(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, mode, taxa_list, taxonomic_level, suffix):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    # keep
    if mode == 'Keep':
        filtered_df = taxon_table_df.loc[taxon_table_df[taxonomic_level].isin(taxa_list)]
    # exclude
    else:
        filtered_df = taxon_table_df.loc[~taxon_table_df[taxonomic_level].isin(taxa_list)]

    ## export table
    export_taxon_table(taxon_table_xlsx, filtered_df, traits_df, metadata_df, suffix)

def sample_filtering(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, mode, samples_list, suffix):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    # keep
    if mode == 'Remove':
        pop_list = [i for i in samples if i in samples_list]
        filtered_df = taxon_table_df.drop(pop_list, axis=1)

    # exclude
    else:
        pop_list = [i for i in samples if i not in samples_list]
        filtered_df = taxon_table_df.drop(pop_list, axis=1)

    ## export table
    export_taxon_table(taxon_table_xlsx, filtered_df, traits_df, metadata_df, suffix)





