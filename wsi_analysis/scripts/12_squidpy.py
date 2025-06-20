import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import seaborn as sns
import numpy as np
import os

path = snakemake.input.anndata_path

sample_id = os.path.basename(path).replace(".h5ad", "")  # Remove extension

label_col = 'labels_0'

adata = sc.read_h5ad(path)

sq.gr.spatial_neighbors(adata, coord_type='generic', radius=140)
sq.gr.nhood_enrichment(adata, cluster_key=label_col)
sq.pl.nhood_enrichment(adata, cluster_key=label_col)

df = pd.DataFrame(adata.uns[f'{label_col}_nhood_enrichment']['zscore'])
df.index = pd.Categorical(adata.obs[label_col].cat.categories)
df.columns = pd.Categorical(adata.obs[label_col].cat.categories)

interaction_df = []
for ct1 in adata.obs[label_col].cat.categories:
    for ct2 in adata.obs[label_col].cat.categories:
        interaction_df.append([sample_id, ct1, ct2, df.loc[ct1, ct2]])
        
interaction_df = pd.DataFrame(interaction_df)
interaction_df.columns = ['sample_id', 'ct1', 'ct2', 'interaction_score']

# Write output
interaction_df.to_csv(snakemake.output.csv_path)
