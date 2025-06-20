import pandas as pd
import numpy as np
import xarray as xr
import spatialproteomics as sp

metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)
files = snakemake.input.zarr_files

# === reading the zarrs into the sample dict ===
sample_dict = {}
sample_to_patient_dict = {}
for i, (sample_id_a, sample_id_b, patient_id) in metadata[["sample_id_a", "sample_id_b", "Histo-Nr"]].iterrows():
    # adding sample a and b to the dictionary if we have images for them
    for sample_id in [sample_id_a, sample_id_b]:
        if pd.isna(sample_id):
            continue
        sample_to_patient_dict[sample_id] = patient_id
        file_mask = [f"{sample_id}." in x for x in files]
        matching_files = [files[i] for i, mask in enumerate(file_mask) if mask]
        if len(matching_files) != 1:
            raise ValueError(f"Matching files for sample {sample_id}: {matching_files}")
        else:
            file = matching_files[0]
            sample_dict[sample_id] = xr.open_zarr(file)
            
patient_to_entity_mapping = dict(zip(metadata['Histo-Nr'].values, metadata['Entity'].values))
entity_grouping_dict = {'LN': 'LN', 'FL 1': 'indolent', 'FL 2': 'indolent', 'FL 3a': 'indolent', 'FL 3b': 'indolent', 'MZL': 'indolent', 'DLBCL': 'aggressive', 'Burkitt': 'aggressive', 'PMBCL': 'aggressive', 'DLBCL': 'aggressive', 'BLBL': 'aggressive'}


# === computing absolute abundances ===
homophily_dfs = []
for sample_id, ds in sample_dict.items():
    nhs = [x for x in ['Neighborhood 3', 'Neighborhood 5'] if x in ds.pp.get_layer_as_df()['_neighborhoods'].unique()]
    homophily_df = pd.DataFrame(ds.nh[nhs].nh.add_neighborhood_obs(features=['degree', 'homophily', 'inter_label_connectivity', 'diversity_index']).pp.get_layer_as_df()[['degree', 'homophily', 'inter_label_connectivity', 'diversity_index']])
    homophily_df['sample_id'] = sample_id
    patient_id = sample_to_patient_dict[sample_id]
    entity = patient_to_entity_mapping[patient_id]
    homophily_df['patient_id'] = patient_id
    homophily_df['entity'] = entity
    homophily_df['entity_class'] = entity_grouping_dict[entity]
    homophily_dfs.append(homophily_df.copy())
    
homophily_df = pd.concat(homophily_dfs, axis=0)

# === storing the results ===
homophily_df.to_csv(snakemake.output.csv_path)
