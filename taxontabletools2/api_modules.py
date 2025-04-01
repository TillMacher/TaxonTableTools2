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

def gbif_accession(taxon_table_xlsx, taxon_table_df, samples, metadata_df, traits_df, higherTaxonKey):
    ## create copies of the dataframes
    taxon_table_df = taxon_table_df.copy()
    metadata_df = metadata_df.copy()
    traits_df = traits_df.copy()

    ## store results here
    species_dict = {'': ['', '', '', '', '']}
    OTU_species = []
    OTU_keys = []
    OTU_genusKeys = []
    OTU_synonyms = []
    OTU_iucnRedListCategory = []
    OTU_habitat = []
    OTU_gbif_link = []

    for OTU in stqdm(taxon_table_df.values.tolist()):
        species = OTU[6]

        ## request GBIF
        if species != '':
            if species not in species_dict.keys():
                query = species.replace(' ', '%20')
                # Initialize a counter for the number of attempts
                attempts = 0
                # Initialize a flag for whether the request was successful
                success = False

                while not success and attempts < 20:
                    try:
                        # Make a GET request to the API
                        response = requests.get(f'https://api.gbif.org/v1/species/match?name={query}')

                        # If the request was successful, parse the response text as JSON
                        if response.status_code == 200:
                            data = json.loads(response.text)
                            success = True
                        else:
                            attempts += 1
                    except requests.exceptions.RequestException as e:
                        # If there was a network problem (e.g. DNS resolution, refused connection, etc), increment the counter and try again
                        attempts += 1

                # If the request was not successful after 20 attempts, print a warning
                if not success:
                    species_dict[species] = [''] * 6
                else:
                    # add data to dict
                    species_dict[species] = [data.get('species', ''), data.get('speciesKey', ''),
                                             data.get('genusKey', ''), data.get('synonym', '')]

                    # Make a GET request to the API to get IUCN Red List Category
                    response = requests.get(
                        f"https://api.gbif.org/v1/species/{data.get('speciesKey', '')}/iucnRedListCategory")
                    if response.status_code == 200:
                        iucn_data = json.loads(response.text)
                        species_dict[species].append(iucn_data.get('code', ''))
                    else:
                        species_dict[species].append('')

                    # Make a GET request to the API to get habitat
                    response = requests.get(
                        f"https://api.gbif.org/v1/species/{data.get('speciesKey', '')}/speciesProfiles")
                    if response.status_code == 200:
                        habitat_data = json.loads(response.text)
                        habitat = ', '.join([j for j in sorted(set([i.get('habitat', '').lower() for i in habitat_data['results']])) if j != ''])
                        species_dict[species].append(habitat)
                    else:
                        species_dict[species].append('')

                time.sleep(0.5)

            # create link
            gbif_url = f'https://www.gbif.org/species/{str(species_dict[species][1])}'

            ## append results
            OTU_species.append(species_dict[species][0])
            OTU_keys.append(str(species_dict[species][1]))
            OTU_genusKeys.append(str(species_dict[species][2]))
            OTU_synonyms.append(str(species_dict[species][3]))
            OTU_iucnRedListCategory.append(species_dict[species][4])
            OTU_habitat.append(species_dict[species][5])
            OTU_gbif_link.append(gbif_url)
        else:
            ## append results
            OTU_species.append('')
            OTU_keys.append('')
            OTU_genusKeys.append('')
            OTU_synonyms.append('')
            OTU_iucnRedListCategory.append('')
            OTU_habitat.append('')
            OTU_gbif_link.append('')

    ## append to traits df
    traits_df['speciesKey'] = OTU_keys
    traits_df['genusKey'] = OTU_genusKeys
    traits_df['GBIF species'] = OTU_species
    traits_df['Synonym'] = OTU_synonyms
    traits_df['iucnRedListCategory'] = OTU_iucnRedListCategory
    traits_df['habitat'] = OTU_habitat
    traits_df['GBIF'] = OTU_gbif_link

    ## export table
    suffix = 'test'
    update_taxon_table(taxon_table_xlsx, taxon_table_df, traits_df, metadata_df, suffix)

