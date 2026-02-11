import xarray as xr
import spatialproteomics as sp
import matplotlib.pyplot as plt
import tifffile
import pandas as pd
import numpy as np
from glob import glob
from tqdm.auto import tqdm
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA

metadata_path = snakemake.input.metadata_path
abundance_df_path = snakemake.input.abundance_df_path
files = snakemake.input.zarr_files
lda_feature_df_path = snakemake.input.lda_feature_df_path

entity_colors = {'LN': '#9dcbec', # light blue
                 'MZL': '#E2AFAC', 
                 'FL 1': '#D59E9C', 'FL 2': '#C98D8B', 'FL 3a': '#BD7C7A', 'FL 3b': '#B16C69', 
                 'DLBCL': '#994A47', 'PMBCL': '#8C3937', 'Burkitt': '#802926', 'BLBL': '#741815'}

entity_class_colors = {'LN': '#9dcbec', 'indolent': '#E2AFAC', 'aggressive': '#741815'}

# red, blue, green, yellow, purple, teal (from the nature color palette: https://www.nature.com/documents/natrev-artworkguide.pdf) 
celltype_colors = {'B cell': '#5799d1', 'T cell': '#ebc850', 'Myeloid cell': '#de6866', 'Dendritic cell': '#4cbcbd', 'Macrophage': '#bb7cb4', 'Stromal cell': '#62b346', 'Endothelial cell': '#bf997d'}

nh_mapping_dict = {'Neighborhood 0': 'T/Dendritic',
                      'Neighborhood 1': 'T',
                      'Neighborhood 2': 'Mixed',
                      'Neighborhood 3': 'B_prol/T_fh',
                      'Neighborhood 4': 'T_h/T_reg',
                      'Neighborhood 5': 'B',
                      'Neighborhood 6': 'Myeloid/Macro'}

# colors for neighborhoods
neighborhood_colors = {'Neighborhood 5': '#5799d1',
                       'Neighborhood 3': '#154e91',
                       'Neighborhood 0': '#fef1c3',
                       'Neighborhood 1': '#ebc850',
                       'Neighborhood 4': '#cca02d',
                       'Neighborhood 6': '#de6866',
                       'Neighborhood 2': '#c8ceda'}

metadata = pd.read_csv(metadata_path, index_col=0)

all_cts = ['B cell', 'B_prol', 'T cell', 'T_h', 'T_h_mem', 'T_h_naive', 'T_reg', 'T_fh', 'T_reg_Helios', 'T_tox', 'T_exh', 'T_progenitor_exh', 'T_terminally_exh', 'T_tox_naive', 'T_tox_mem', 'Myeloid cell', 'Dendritic cell', 'Macrophage', 'M2', 'Stromal cell', 'Endothelial cell']
all_nhs = list(nh_mapping_dict.values())

entity_grouping_dict = {'LN': 'LN', 'FL 1': 'indolent', 'FL 2': 'indolent', 'FL 3a': 'indolent', 'FL 3b': 'indolent', 'MZL': 'indolent', 'DLBCL': 'aggressive', 'Burkitt': 'aggressive', 'PMBCL': 'aggressive', 'DLBCL': 'aggressive', 'BLBL': 'aggressive'}

# for this figure, I want to sort by entity and B cell content, but also keep replicates next to one another
# adding the celltype abundances which have been precomputed in the snakemake pipeline
abundance_df = pd.read_csv(abundance_df_path, index_col=0)
# normalizing the cell type abundances per sample
# normalizing per sample so that each sample sums up to one
abundance_df = abundance_df.div(abundance_df.sum(axis=1), axis=0)
# merging with the metadata
df = pd.merge(metadata, abundance_df, left_on='Histo-Nr', right_index=True)

# Convert 'entity' column to a categorical type with the custom order
df['Entity'] = pd.Categorical(df['Entity'], categories=entity_colors.keys(), ordered=True)
df = df.sort_values(by=['Entity', 'B cell'])

# === reading the zarrs into the sample dict ===
sample_dict = {}
sample_to_patient_dict = {}
for i, (sample_id_a, sample_id_b, patient_id) in tqdm(df[["sample_id_a", "sample_id_b", "Histo-Nr"]].iterrows()):
    # adding sample a and b to the dictionary if we have images for them
    for sample_id in [sample_id_a, sample_id_b]:
        if pd.isna(sample_id):
            continue
        sample_to_patient_dict[sample_id] = patient_id
        file_mask = [f"{sample_id}.zarr" in x for x in files]
        # print(files)
        # file_mask = [sample_id == x for x in files]
        matching_files = [files[i] for i, mask in enumerate(file_mask) if mask]
        if len(matching_files) != 1:
            raise ValueError(f"Matching files for sample {sample_id}: {matching_files}")
        else:
            file = matching_files[0]
            ds = xr.open_zarr(file)
            # also renaming the neighborhoods here
            nh_mapping_dict_tmp  = {k: v for k, v in nh_mapping_dict.items() if k in ds.pp.get_layer_as_df()['_neighborhoods'].unique()}
            sample_dict[sample_id] = ds.nh.set_neighborhood_colors(neighborhood_colors.keys(), neighborhood_colors.values()).nh.set_neighborhood_name(list(nh_mapping_dict_tmp.keys()), list(nh_mapping_dict_tmp.values()))

                        
patient_to_entity_mapping = dict(zip(df['Histo-Nr'].values, df['Entity'].values))

# === computing the cell sizes df ===
dfs = []
for sample, ds in tqdm(sample_dict.items()):
    df_tmp = ds.pp.add_observations('area').pp.get_layer_as_df()[['area', 'degree']]
    df_tmp['num_cells'] = df_tmp.shape[0]
    df_tmp['sample'] = sample
    df_tmp['entity'] = patient_to_entity_mapping[sample_to_patient_dict[sample]]
    dfs.append(df_tmp)
df_degree_size_tmp = pd.concat(dfs)
                                                 
# Aggregating the number of cells and area by entity
df_degree_size = df_degree_size_tmp.groupby("entity").agg({
    "area": "mean",       # Mean of areas
    "num_cells": "mean"
}).reset_index()

# Add the number of cells directly in the aggregation
# df_degree_size["num_cells"] = df_degree_size_tmp.groupby("entity").size().values

# Sort entities
df_degree_size["entity"] = pd.Categorical(
    df_degree_size["entity"], 
    categories=entity_colors.keys(), 
    ordered=True
)
df_degree_size = df_degree_size.sort_values("entity").reset_index(drop=True)
df_degree_size['entity_class'] = [entity_grouping_dict[x] for x in df_degree_size['entity']]

# === LDA ===
# constructing the df for the LDA
# adding neighborhood abundances and network features (could also consider adding things like number of cells, mean cell area, etc.)
# this is done per sample and not per patient, hence we reconstruct the ct abundances here and don't use the ones computed by the snakemake pipeline

dfs = []
for sample, ds in tqdm(sample_dict.items()):
    df_tmp = []
    num_cells = ds.sizes['cells']
    
    # step 1: get the ct abundances
    ct_abundances = ds.pp.get_layer_as_df('_la_layers')['labels_3'].value_counts()
    
    for ct in all_cts:
        df_tmp.append(ct_abundances.get(ct, 0) / num_cells)
        
    # step 2: get the neighborhood abundances
    nh_abundances = ds.pp.get_layer_as_df('_obs')['_neighborhoods'].value_counts()
    for nh in all_nhs:
        df_tmp.append(nh_abundances.get(nh, 0) / num_cells)
        
    # step 3: get network-based features
    # the network level features are probably not good enough, will instead compute them per cell and then aggregate by taking the mean
    # these are already computed, no need to recompute
    network_based_features = ['degree', 'homophily', 'inter_label_connectivity', 'diversity_index']
    network_features = list(ds.pp.get_layer_as_df('_obs')[network_based_features].mean(axis=0).values)
    df_tmp.extend(network_features)
    
    # step 4: get morphology-based features (area, number of cells)
    # number of cells
    df_tmp.append(num_cells)
    # eccentricity is 0 if the ellipse is a 0 and goes up to 1
    df_tmp.extend(list(ds.pp.add_observations(['area', 'eccentricity']).pp.get_layer_as_df('_obs')[['area', 'eccentricity']].mean(axis=0).values))
    
    # adding the sample and entity
    patient_id = sample_to_patient_dict[sample]
    entity = patient_to_entity_mapping[patient_id]
    df_tmp.extend([sample, patient_id, entity])
    dfs.append(df_tmp)
    
df_final = pd.DataFrame(dfs)
df_final.columns = all_cts + all_nhs + ['Degree (mean)', 'Homophily (mean)', 'Inter-label connectivity (mean)', 'Diversity index (mean)', 'Number of cells', 'Cell area (mean)', 'Cell eccentricity (mean)', 'sample_id', 'patient_id', 'entity']
df_final['entity_class'] = [entity_grouping_dict[x] for x in df_final['entity']]

df_final.to_csv(lda_feature_df_path)

# running LDA
# Specify columns to exclude from scaling
exclude_cols = ['sample_id', 'patient_id', 'entity', 'entity_class']
scale_cols = [col for col in df_final.columns if col not in exclude_cols]

# Center and scale the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df_final[scale_cols])

# Encode the target variable

label_encoder = LabelEncoder()
entity_labels = label_encoder.fit_transform(df_final['entity_class'])

# Run LDA
lda = LDA(n_components=2)  # Use 2 components for 2D plotting
lda_components = lda.fit_transform(scaled_data, entity_labels)
lda_df = pd.DataFrame(lda_components, columns=['LD1', 'LD2'])
lda_df['entity_class'] = df_final['entity_class']

# Explained variance ratio for LDA
explained_variance = lda.explained_variance_ratio_ * 100

# Prepare loadings and sort them
# Prepare loadings dynamically based on the number of components
loadings = pd.DataFrame(lda.scalings_, 
                        columns=[f'LD{i+1}' for i in range(lda.scalings_.shape[1])], 
                        index=scale_cols)

# Sort loadings for LD1 and LD2
loadings = loadings.sort_values(by='LD1', ascending=False)

# === Neighborhood composition ===
# computing the neighborhood composition
image_container = sp.ImageContainer(sample_dict)
nh_composition_df = image_container.get_neighborhood_composition()
df = df.reset_index(drop=True)

# === exporting the data ===
nh_composition_df.to_csv(snakemake.output.neighborhood_composition_path)
df.to_csv(snakemake.output.celltype_abundance_df)
loadings.to_csv(snakemake.output.lda_loadings_path)
lda_df.to_csv(snakemake.output.lda_path)
df_degree_size.to_csv(snakemake.output.degree_size_df_path)
pd.DataFrame(explained_variance).to_csv(snakemake.output.lda_ev_path)