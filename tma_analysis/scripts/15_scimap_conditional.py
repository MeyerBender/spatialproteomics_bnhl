import scimap as sm
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

# Read input from Snakemake
adata = sc.read_h5ad(snakemake.input.anndata_path)
adata.obs_names_make_unique()
label_col = snakemake.wildcards.label  # Get label from wildcard

# === Spatial interaction ===
adata = sm.tl.spatial_interaction(
    adata, 
    x_coordinate='centroid-1', y_coordinate='centroid-0', phenotype=label_col,
    method='delaunay', 
    permutation=100,
    imageid='sample_id',
    pval_method='zscore',
    normalization='conditional',
    label='scimap_delaunay_conditional_zscore'
)

# Write output
adata.write(snakemake.output.anndata_path)
