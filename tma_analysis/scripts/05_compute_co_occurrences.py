import xarray as xr
import spatialproteomics as sp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from glob import glob
from tqdm.auto import tqdm
import squidpy as sq

ds = xr.open_zarr(snakemake.input.zarr_path)
metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)
sample_id = snakemake.input.zarr_path.split('/')[-1].replace('.zarr', '')
sample_id_list = sample_id.split('_')
sample_id = f'{sample_id_list[0]}_{sample_id_list[1]}_{sample_id_list[2]}'

abundance_df = pd.read_csv("/g/huber/users/meyerben/data/codex/BNHL/ct_abundances/ct_abundances_level_0.csv", index_col=0)
# normalizing the cell type abundances per sample
# normalizing per sample so that each sample sums up to one
abundance_df = abundance_df.div(abundance_df.sum(axis=1), axis=0)
# merging with the metadata
df = pd.merge(metadata, abundance_df, left_on='Histo-Nr', right_index=True)
patient_to_entity_mapping = dict(zip(df['Histo-Nr'].values, df['Entity'].values))

entity_colors = {'LN': '#9dcbec', # light blue
                 'CLL': '#FAD0CE', 'LPL': '#EEBFBD', 'MZL': '#E2AFAC', 
                 'FL 1': '#D59E9C', 'FL 2': '#C98D8B', 'FL 3a': '#BD7C7A', 'FL 3b': '#B16C69', 
                 'MCL': '#A55B58', 'DLBCL': '#994A47', 'PMBCL': '#8C3937', 'Burkitt': '#802926', 'BLBL': '#741815'}

# Convert 'entity' column to a categorical type with the custom order
df['Entity'] = pd.Categorical(df['Entity'], categories=entity_colors.keys(), ordered=True)
df = df.sort_values(by=['Entity', 'B cell'])

# === reading the zarrs into the sample dict ===
sample_dict = {}
sample_to_patient_dict = {}
for i, (sample_id_a, sample_id_b, patient_id) in tqdm(df[["sample_id_a", "sample_id_b", "Histo-Nr"]].iterrows()):
    # adding sample a and b to the dictionary if we have images for them
    for sample_id_tmp in [sample_id_a, sample_id_b]:
        if pd.isna(sample_id_tmp):
            continue
        sample_to_patient_dict[sample_id_tmp] = patient_id


all_cts = ['B cell', 'B_prol', 'T_h', 'T_h_mem', 'T_h_naive', 'T_reg', 'T_fh', 'T_reg_Helios', 'T_tox', 'T_exh', 'T_progenitor_exh', 'T_terminally_exh', 'T_tox_naive', 'T_tox_mem', 'Myeloid cell', 'Dendritic cell', 'Macrophage', 'M2', 'Stromal cell', 'Endothelial cell']

def get_probability_df_per_sample(sample_id, ct):
    # for some samples that did not pass QC, we cannot to return anything, because we don't have expression matrices for those
    if sample_id not in sample_to_patient_dict.keys():
        return None
    try:
        adata = ds.tl.convert_to_anndata()
    except AssertionError:
        return None
    adata.obs['_labels'] = adata.obs['_labels'].astype('category')

    # computing co-occurrence (from 5 to 150 microns)
    sq.gr.co_occurrence(adata, cluster_key="_labels", interval=np.array(np.arange(0, 510, 10)))
    
    occurrence_data = adata.uns['_labels_co_occurrence']
    out = occurrence_data["occ"]
    interval = occurrence_data["interval"][1:]
    categories = adata.obs['_labels'].cat.categories

    # in case the cell type does not appear in the sample, we do not add anything to the dataframe
    try:
        idx = np.where(categories == ct)[0][0]
        df = pd.DataFrame(out[idx, :, :].T, columns=categories).melt(var_name='_labels', value_name="probability")
        df["distance"] = np.tile(interval, len(categories))
        df.columns = ['celltype', 'probability', 'distance']

        # adding sample ID, patient ID, and entity to the dataframe
        patient_id = sample_to_patient_dict[sample_id]
        df['sample_id'] = sample_id
        df['patient_id'] = patient_id
        df['entity'] = patient_to_entity_mapping[patient_id]

        # also adding the conditional in there
        df["conditional_celltype"] = ct

        return df
    except IndexError:
        return None

dfs = []
for ct in all_cts:
    probability_df = get_probability_df_per_sample(sample_id=sample_id, ct=ct)
    if probability_df is not None:
        dfs.append(probability_df)

if dfs:  # Check if the list is not empty
    df = pd.concat(dfs, axis=0).reset_index(drop=True)
else:  # Create an empty DataFrame
    df = pd.DataFrame()

df.to_csv(snakemake.output.co_occurrence_df)
