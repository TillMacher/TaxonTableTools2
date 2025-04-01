import math
import time
import pandas as pd
import numpy as np
from pathlib import Path
from taxontabletools2.utilities import add_traits, filter_taxontable, export_taxon_table, update_taxon_table
import plotly.graph_objects as go
import plotly.express as px
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
import streamlit as st
import plotly.colors
from taxontabletools2.utilities import load_df, collect_traits, strip_traits
import os, subprocess
from playwright.sync_api import sync_playwright

def convert_to_perlodes(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, tool_settings):
    # Make copies of input dataframes to prevent altering the originals
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    river_types_dict = {i:j for i,j in metadata_df[['Sample', 'Perlodes_river_type']].values.tolist()}
    taxa_list_dict = {i:j for i,j in metadata_df[['Sample', 'Perlodes_taxa_list']].values.tolist()}
    usage_dict = {i:j for i,j in metadata_df[['Sample', 'Perlodes_usage']].values.tolist()}

    # Ensure all samples have metadata in all three categories
    samples = [
        i for i in samples
        if i in river_types_dict and pd.notna(river_types_dict[i]) and
           i in taxa_list_dict and pd.notna(taxa_list_dict[i]) and
           i in usage_dict and pd.notna(usage_dict[i])
            ]

    if len(samples) == 0:
        st.warning('Please fill out the required metadata for at least one sample!')
        return

    # Extract relevant settings from the tool settings
    presence_absence = tool_settings['presence_absence']
    taxon_table_taxonomy = taxon_table_df.columns.tolist()[1:7]

    all_reads = []
    for sample in samples:
        total_reads = taxon_table_df[sample].sum()
        taxon_dict = {}
        for OTU in taxon_table_df[taxon_table_taxonomy + [sample]].values.tolist():
            taxon = [item for item in OTU[0:6] if item != ''][-1]
            if taxon not in taxon_dict.keys():
                n_reads = OTU[-1]
                taxon_dict[taxon] = n_reads
            else:
                n_reads = OTU[-1]
                taxon_dict[taxon] = taxon_dict[taxon] + n_reads
        if presence_absence == True:
            pa_list = [1 if i != 0 else 0 for i in taxon_dict.values()]
            all_reads.append(pa_list)
        elif presence_absence == 'Relative':
            rel_list = [round(i / total_reads * 100, 2) for i in taxon_dict.values()]
            all_reads.append(rel_list)
        else:
            all_reads.append(list(taxon_dict.values()))

    # create initial df
    index = [i for i in range(1,1+len(taxon_dict.keys()))]
    species = list(taxon_dict.keys())
    perlodes_df = pd.DataFrame(species, columns=['species'])
    perlodes_df.insert(0, '', index)
    perlodes_df = pd.concat([perlodes_df, pd.DataFrame(all_reads, index=samples).transpose()], axis=1)

    # convert the df to match the (complicated) perlodes input format
    columns = ['ID_ART', 'TAXON_NAME'] + samples
    gewässertyp = ['Gewässertyp', ''] + [river_types_dict[i] for i in samples]
    Taxaliste = ['Taxaliste', ''] + [taxa_list_dict[i] for i in samples]
    Nutzung = ['Nutzung', ''] + [usage_dict[i] for i in samples]
    taxa = [[1] + i[1:] for i in perlodes_df.values.tolist()]
    df_list = [gewässertyp] + [Taxaliste] + [Nutzung] + taxa
    perlodes_input_df = pd.DataFrame(df_list, columns=columns)

    # write the filtered list to a dataframe
    if presence_absence == True:
        perlodes_directory = Path(str(path_to_outdirs) + "/" + "Perlodes" + "/" + taxon_table_xlsx.stem)
        perlodes_xlsx = Path(str(perlodes_directory) + "_Perlodes_PA.xlsx")
        perlodes_input_df.to_excel(perlodes_xlsx, index=False)
    elif presence_absence == 'Relative':
        perlodes_directory = Path(str(path_to_outdirs) + "/" + "Perlodes" + "/" + taxon_table_xlsx.stem)
        perlodes_xlsx = Path(str(perlodes_directory) + "_Perlodes_RELATIVE_ABUNDANCE.xlsx")
        perlodes_input_df.to_excel(perlodes_xlsx, index=False)
    else:
        perlodes_directory = Path(str(path_to_outdirs) + "/" + "Perlodes" + "/" + taxon_table_xlsx.stem)
        perlodes_xlsx = Path(str(perlodes_directory) + "_Perlodes_ABUNDANCE.xlsx")
        perlodes_input_df.to_excel(perlodes_xlsx, index=False)

    st.success(f'Wrote Perlodes input file to: {perlodes_xlsx}')

def convert_to_phylib(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, tool_settings):

    #get the taxonomy from the operational taxon list
    operational_taxon_list_df = pd.read_excel(Path(operational_taxon_list), sheet_name="TTT import").fillna('')

    # load the taxon table and create a list
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    meta_data_df = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)
    TaXon_table_taxonomy = TaXon_table_df.columns.tolist()[0:7]
    samples_list = TaXon_table_df.columns.tolist()[10:]

    ## load the metadata -> freshwater type
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
    metadata_loc = Meta_data_table_df.columns.tolist().index(meta_data_to_test)
    types_dict = {i[0]:i[1] for i in Meta_data_table_df[['Samples', meta_data_to_test]].values.tolist()}

    ## drop samples with metadata called nan (= empty)
    drop_samples = [i[0] for i in Meta_data_table_df.values.tolist() if i[metadata_loc] == ""]

    ## test if samples have metadata
    if drop_samples != []:
        sg.PopupError("Please fill out all the metadata for all samples.")

    ## test if all metadata for Phylib is available
    elif len(set([True if i in Meta_data_table_df.columns.tolist() else False for i in ['Ökoregion', 'Makrophytenverödung', 'Begründung', 'Helophytendominanz', 'Diatomeentyp', 'Phytobenthostyp', 'Makrophytentyp', 'WRRL-Typ', 'Gesamtdeckungsgrad']])) != 1:
        sg.PopupError("Please fill out all the required phylib metadata for all samples.")

    elif sorted(Meta_data_table_samples) == sorted(samples_list):
        # store hits and dropped OTUs
        hit_list, dropped_list = [], []

        # loop through the taxon table
        for taxonomy in TaXon_table_df[TaXon_table_taxonomy].values.tolist():
            ## test species
            if taxonomy[6] != '' and taxonomy[6] in operational_taxon_list_df['Species'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Species'].str.contains(taxonomy[6])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test genus
            elif taxonomy[5] != '' and taxonomy[5] in operational_taxon_list_df['Genus'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Genus'].str.contains(taxonomy[5])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test family
            elif taxonomy[4] != '' and taxonomy[4] in operational_taxon_list_df['Family'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Family'].str.contains(taxonomy[4])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test order
            elif taxonomy[3] != '' and taxonomy[3] in operational_taxon_list_df['Order'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Order'].str.contains(taxonomy[3])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test class
            elif taxonomy[2] != '' and taxonomy[2] in operational_taxon_list_df['Class'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Class'].str.contains(taxonomy[2])][['DV-NR.', 'Taxon']].values.tolist()[0]

            ## test phylum
            elif taxonomy[1] != '' and taxonomy[1] in operational_taxon_list_df['Phylum'].values.tolist():
                ## merge the OTU's taxonomy with OTL information
                res = taxonomy + operational_taxon_list_df[operational_taxon_list_df['Phylum'].str.contains(taxonomy[1])][['DV-NR.', 'Taxon']].values.tolist()[0]
            else:
                res = taxonomy + ['', '']
                dropped_list.append(taxonomy)

            hit_list.append(res)

        ## create a dataframe
        ## create a new df and export it as TaXon table
        df = pd.DataFrame(hit_list, columns=TaXon_table_df.columns[0:7].values.tolist() + ['DV-NR.', 'Taxon (Phylib)'])
        concatenated_df = pd.concat([df, TaXon_table_df[samples_list + ['Similarity', 'Status', 'seq']]], axis=1)
        reordered_df = concatenated_df[['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'Status', 'DV-NR.', 'Taxon (Phylib)', 'seq'] + samples_list]
        reordered_metadata_df = add_metadata(reordered_df, meta_data_df)
        phylib_TaXon_Table = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem + "_phylib.xlsx")
        reordered_metadata_df.to_excel(phylib_TaXon_Table, index=False)

        # create an output list for phylib
        messwerte_list = []
        messtellen_list = []

        for sample in samples_list:
            ## collect information about the sampling site
            sample_metadata_df = Meta_data_table_df.loc[Meta_data_table_df['Samples'] == sample]
            messstelle = sample
            ökoregion = sample_metadata_df['Ökoregion'].values.tolist()[0]
            Makrophytenverödung = sample_metadata_df['Makrophytenverödung'].values.tolist()[0]
            Begründung = sample_metadata_df['Begründung'].values.tolist()[0]
            Helophytendominanz = sample_metadata_df['Helophytendominanz'].values.tolist()[0]
            Diatomeentyp = sample_metadata_df['Diatomeentyp'].values.tolist()[0]
            Phytobenthostyp = sample_metadata_df['Phytobenthostyp'].values.tolist()[0]
            Makrophytentyp = sample_metadata_df['Makrophytentyp'].values.tolist()[0]
            WRRL_Typ = sample_metadata_df['WRRL-Typ'].values.tolist()[0]
            Gesamtdeckungsgrad = sample_metadata_df['Gesamtdeckungsgrad'].values.tolist()[0]
            messtellen_list.append([messstelle, ökoregion, Makrophytenverödung, Begründung, Helophytendominanz, Diatomeentyp, Phytobenthostyp, Makrophytentyp, WRRL_Typ, Gesamtdeckungsgrad])

            ## calcualte sum of reads for taxa with multiple OTUs
            sample_df = reordered_metadata_df[[sample, 'DV-NR.', 'Taxon (Phylib)']]
            reads_dict = {}
            for taxon in sample_df.values.tolist():
                key = taxon[2]
                values = taxon[:2]
                if key != '':
                    if key not in reads_dict.keys():
                        reads_dict[key] = values
                    else:
                        reads_dict[key] = [reads_dict[key][0] + values[0], values[1]]

            ## remove duplicates
            samples_taxa = [[key]+values for key,values in reads_dict.items() if values[0] != 0]

            ## if PA data is required: Convert to 1/0
            if presence_absence == True:
                samples_taxa = [[i[0], 1, i[2]] for i in samples_taxa]

            ## calculate the overall number of reads/specimens. This is required for later relative abundance calculation
            sum_measurement = sum([i[1] for i in samples_taxa])

            ## loop through all taxa
            for taxon in samples_taxa:
                ## if the taxon (assigned to the OTU) is present in the sample and present on the OTL, continue
                ## add all relevant information to list (for df)
                probe = sample
                taxon_id = taxon[2]
                taxonname = taxon[0]
                form = "o.A."
                messwert = taxon[1] / sum_measurement * 100
                einheit = "%"
                cf = ""
                messwerte_list.append([messstelle, probe, taxon_id, taxonname, form, messwert, einheit, cf])

        phylib_df_1 = pd.DataFrame(messtellen_list, columns=["Messstelle", "Ökoregion", "Makrophytenverödung", "Begründung", "Helophytendominanz", "Diatomeentyp", "Phytobenthostyp", "Makrophytentyp", "WRRL-Typ", "Gesamtdeckungsgrad"])
        phylib_df_2 = pd.DataFrame(messwerte_list, columns=["Messstelle", "Probe", "Taxon", "Taxonname", "Form", "Messwert", "Einheit", "cf"])

        if presence_absence == False:
            phylib_directory = Path(str(path_to_outdirs) + "/" + "Phylib" + "/" + TaXon_table_xlsx.stem)
            phylib_xlsx = Path(str(phylib_directory) + "_phylib_ABUNDANCE.xlsx")
            writer = pd.ExcelWriter(phylib_xlsx, engine='xlsxwriter')
        else:
            phylib_directory = Path(str(path_to_outdirs) + "/" + "Phylib" + "/" + TaXon_table_xlsx.stem)
            phylib_xlsx = Path(str(phylib_directory) + "_phylib_PA.xlsx")
            writer = pd.ExcelWriter(phylib_xlsx, engine='xlsxwriter')
        phylib_df_1.to_excel(writer, sheet_name='Messstelle', index=False)
        phylib_df_2.to_excel(writer, sheet_name='Messwerte', index=False)
        writer.save()

        ################################################################################################################
        ## create some plots and provide statistics
        unique_original_taxa = list(k for k, _ in itertools.groupby(reordered_metadata_df[['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].values.tolist()))
        len(unique_original_taxa)

        unique_OTL_taxa = [i for i in set(reordered_metadata_df['Taxon (Phylib)'].values.tolist()) if i != '']
        len(unique_OTL_taxa)

        phylib_conversion_loss = {}
        for taxon in unique_OTL_taxa:
            all_taxa = reordered_metadata_df[reordered_metadata_df['Taxon (Phylib)'].str.contains(taxon)][['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].values.tolist()
            unique_taxa = set([''.join(i) for i in all_taxa])
            phylib_conversion_loss[taxon] = len(unique_taxa)

        perlodes_conversion_loss = {k: v for k, v in sorted(phylib_conversion_loss.items(), key=lambda item: item[1])}

        fig = go.Figure()
        x_values = [i[:19] for i in list(perlodes_conversion_loss.keys())[-30:][::-1]] ## [:20] to cut names off
        y_values = list(perlodes_conversion_loss.values())[-30:][::-1]
        fig.add_trace(go.Bar(x=x_values, y=y_values, marker_color='blue'))
        fig.update_yaxes(title='# taxa')
        fig.update_layout(title='', showlegend=False, font_size=14, width=width_value, height=height_value, template=template)
        phylib_plot = Path(str(phylib_directory) + "_phylib.pdf")
        fig.write_image(phylib_plot)
        phylib_plot = Path(str(phylib_directory) + "_phylib.html")
        fig.write_html(phylib_plot)

def convert_to_diathor(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, tool_settings):
    # Make copies of input dataframes to prevent altering the originals
    taxon_table_df = taxon_table_df.copy()

    # Extract relevant settings from the tool settings
    presence_absence = tool_settings['presence_absence']
    taxon_table_taxonomy = taxon_table_df.columns.tolist()[1:7]

    all_reads = []
    for sample in samples:
        total_reads = taxon_table_df[sample].sum()
        taxon_dict = {}
        for OTU in taxon_table_df[taxon_table_taxonomy + [sample]].values.tolist():
            taxon = [item for item in OTU[0:6] if item != ''][-1]
            if taxon not in taxon_dict.keys():
                n_reads = OTU[-1]
                taxon_dict[taxon] = n_reads
            else:
                n_reads = OTU[-1]
                taxon_dict[taxon] = taxon_dict[taxon] + n_reads
        if presence_absence == True:
            pa_list = [1 if i != 0 else 0 for i in taxon_dict.values()]
            all_reads.append(pa_list)
        elif presence_absence == 'Relative':
            rel_list = [round(i / total_reads * 100, 2) for i in taxon_dict.values()]
            all_reads.append(rel_list)
        else:
            all_reads.append(list(taxon_dict.values()))

    index = [i for i in range(1,1+len(taxon_dict.keys()))]
    species = list(taxon_dict.keys())

    diathor_df = pd.DataFrame(species, columns=['species'])
    diathor_df.insert(0, '', index)
    diathor_df = pd.concat([diathor_df, pd.DataFrame(all_reads, index=samples).transpose()], axis=1)

    # write the filtered list to a dataframe
    if presence_absence == True:
        diathor_directory = Path(str(path_to_outdirs) + "/" + "Diathor" + "/" + taxon_table_xlsx.stem)
        diathor_csv = Path(str(diathor_directory) + "_DiaThor_PA.csv")
        diathor_df.to_csv(diathor_csv, index=False, sep=',')
    elif presence_absence == 'Relative':
        diathor_directory = Path(str(path_to_outdirs) + "/" + "Diathor" + "/" + taxon_table_xlsx.stem)
        diathor_csv = Path(str(diathor_directory) + "_DiaThor_RELATIVE_ABUNDANCE.csv")
        diathor_df.to_csv(diathor_csv, index=False, sep=',')
    else:
        diathor_directory = Path(str(path_to_outdirs) + "/" + "DiaThor" + "/" + taxon_table_xlsx.stem)
        diathor_csv = Path(str(diathor_directory) + "_DiaThor_ABUNDANCE.csv")
        diathor_df.to_csv(diathor_csv, index=False, sep=',')

    st.success(f'Wrote Diathor input file to: {diathor_csv}')

    ## Calculate Diathor
    output_folder = diathor_csv.parent.joinpath(diathor_csv.stem)
    os.makedirs(output_folder, exist_ok=True)
    output_xlsx = output_folder.joinpath(str(diathor_csv.name).replace('.csv', '.xlsx'))
    log_txt = Path(str(output_xlsx).replace('.xlsx', '.log'))
    f = open(log_txt, 'w')
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = Path(os.path.join(script_dir, 'Rscripts', 'Diathor.R'))
    subprocess.call(['Rscript', script_path, '-i', diathor_csv, '-o', output_xlsx], stdout=f)
    f.close()
    st.success(f'Wrote Diathor calculations to: {output_xlsx}')

def convert_to_efi(path_to_outdirs, taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, tool_settings):
    # Make copies of input dataframes to prevent altering the originals
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()

    # Extract relevant settings from the tool settings
    presence_absence = tool_settings['presence_absence']
    taxon_table_taxonomy = taxon_table_df.columns.tolist()[1:7]

    # Prepare metadata
    efi_metadata_df = metadata_df[['Sample'] + tool_settings['efi_columns']]

    # Filter samples
    samples = [sample for sample in samples if '' not in efi_metadata_df[efi_metadata_df['Sample'] == sample].values.tolist()[0]]

    if len(samples) == 0:
        st.warning('Please fill out the required metadata for at least one sample!')
        return

    # collect species
    i = 1
    efi_df_values = []
    for sample in samples:
        if presence_absence == 'Genetic Diversity':
            taxa = [i for i in taxon_table_df[taxon_table_df[sample] != 0]['Species'].values.tolist() if i != '']
            taxon_dict = {i:taxa.count(i) for i in set(taxa)}
        else:
            taxon_dict = {}
            total_reads = taxon_table_df[sample].sum()
            for OTU in taxon_table_df[taxon_table_taxonomy + [sample]].values.tolist():
                taxon = [item for item in OTU[0:6] if item != ''][-1]
                if taxon not in taxon_dict.keys():
                    n_reads = OTU[-1]
                    taxon_dict[taxon] = n_reads
                else:
                    n_reads = OTU[-1]
                    taxon_dict[taxon] = taxon_dict[taxon] + n_reads

        efi_sample_metadata = efi_metadata_df[efi_metadata_df['Sample'] == sample].values.tolist()
        for taxon, reads in taxon_dict.items():
            if reads != 0:
                if presence_absence == True:
                    total_number_run1 = 2
                    number_length_below_150 = 1
                    number_length_above_150 = 1
                elif presence_absence == 'Relative':
                    total_number_run1 = math.ceil(reads / total_reads * 100)
                    number_length_below_150 = math.ceil(total_number_run1 / 2)
                    number_length_above_150 = math.ceil(total_number_run1 / 2)
                    total_number_run1 = number_length_above_150 * 2 # needs to be updated to fit the sum of both
                else:
                    number_length_below_150 = math.ceil(reads/2)
                    number_length_above_150 = math.ceil(reads/2)
                    total_number_run1 = number_length_above_150 * 2 # needs to be updated to fit the sum of both
                species_values = efi_sample_metadata[0] + [0, taxon, total_number_run1, number_length_below_150, number_length_above_150, 'abc', i, taxon]
                efi_df_values.append(species_values)
                i +=1

    # construct EFI dataframe
    efi_df = pd.DataFrame(efi_df_values, columns=['Sample.code'] + [i.replace('EFI_', '') for i in tool_settings['efi_columns']] + ['Medit', 'Species', 'Total.number.run1', 'Number.length.below.150', 'Number.length.over.150', 'Sampling.location', 'code', 'species'])

    # add dummy values
    # other EFI will crash in some cases...
    efi_table = Path(os.path.join(st.session_state['script_dir'], 'WFD_conversion', 'efi_table.xlsx'))
    efi_dummy_df = pd.read_excel(efi_table, sheet_name='Sheet3').fillna('')
    efi_dummy_df = efi_dummy_df[efi_df.columns.tolist()]
    efi_df = pd.concat([efi_df, efi_dummy_df], ignore_index=True)
    efi_df['code'] = [i for i in range(1, len(efi_df)+1)]

    # write the filtered list to a dataframe
    if presence_absence == True:
        efi_directory = Path(str(path_to_outdirs) + "/" + "EFI" + "/" + taxon_table_xlsx.stem)
        efi_xlsx = Path(str(efi_directory) + "_EFI_PA.xlsx")
        efi_df.to_excel(efi_xlsx, index=False)
    elif presence_absence == 'Relative':
        efi_directory = Path(str(path_to_outdirs) + "/" + "EFI" + "/" + taxon_table_xlsx.stem)
        efi_xlsx = Path(str(efi_directory) + "_EFI_RELATIVE_ABUNDANCE.xlsx")
        efi_df.to_excel(efi_xlsx, index=False)
    elif presence_absence == 'Genetic Diversity':
        efi_directory = Path(str(path_to_outdirs) + "/" + "EFI" + "/" + taxon_table_xlsx.stem)
        efi_xlsx = Path(str(efi_directory) + "_GENETIC_DIVERSITY.xlsx")
        efi_df.to_excel(efi_xlsx, index=False)
    else:
        efi_directory = Path(str(path_to_outdirs) + "/" + "EFI" + "/" + taxon_table_xlsx.stem)
        efi_xlsx = Path(str(efi_directory) + "_EFI_ABUNDANCE.xlsx")
        efi_df.to_excel(efi_xlsx, index=False)

    st.success(f'Wrote EFI input file to: {efi_xlsx}')

    ## Calculate EFI
    output_folder = efi_xlsx.parent.joinpath(efi_xlsx.stem)
    os.makedirs(output_folder, exist_ok=True)
    output_xlsx = output_folder.joinpath(str(efi_xlsx.name).replace('.xlsx', '_EFI.xlsx'))
    log_txt = Path(str(output_xlsx).replace('.xlsx', '.log'))
    f = open(log_txt, 'w')
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = Path(os.path.join(script_dir, 'Rscripts', 'EFI.R'))
    subprocess.call(['Rscript', script_path, '-i', efi_xlsx, '-o', output_xlsx], stdout=f)
    f.close()
    st.success(f'Wrote EFI calculations to: {output_xlsx}')

