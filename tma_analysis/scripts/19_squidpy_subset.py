import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import seaborn as sns
import numpy as np
import os

path = snakemake.input.anndata_path
metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)

sample_id = os.path.basename(path).replace(".h5ad", "")  # Remove extension

# Find matching metadata row
matching_rows = metadata[(metadata["sample_id_a"] == sample_id) | (metadata["sample_id_b"] == sample_id)]

# Ensure exactly one match is found
if len(matching_rows) != 1:
    raise ValueError(f"Expected exactly one metadata row for sample_id '{sample_id}', but found {len(matching_rows)}")

# Get metadata as a dictionary (excluding sample_id_a and sample_id_b)
metadata_dict = matching_rows.drop(columns=["sample_id_a", "sample_id_b"]).iloc[0].to_dict()

label_col = snakemake.wildcards.label  # Get label from wildcard

adata = sc.read_h5ad(path)
adata = adata[adata.obs['_neighborhoods'].isin(['Neighborhood 3', 'Neighborhood 5'])].copy()
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=140)
sq.gr.nhood_enrichment(adata, cluster_key=label_col)
sq.pl.nhood_enrichment(adata, cluster_key=label_col)

df = pd.DataFrame(adata.uns[f'{label_col}_nhood_enrichment']['zscore'])
df.index = pd.Categorical(adata.obs[label_col].cat.categories)
df.columns = pd.Categorical(adata.obs[label_col].cat.categories)

interaction_df = []
for ct1 in adata.obs[label_col].cat.categories:
    for ct2 in adata.obs[label_col].cat.categories:
        interaction_df.append([sample_id, metadata_dict['Histo-Nr'], metadata_dict['Entity'], metadata_dict['Entity_Class'], ct1, ct2, df.loc[ct1, ct2]])
        
interaction_df = pd.DataFrame(interaction_df)
interaction_df.columns = ['sample_id', 'patient', 'entity', 'entity_class', 'ct1', 'ct2', 'interaction_score']

# Write output
interaction_df.to_csv(snakemake.output.csv_path)
