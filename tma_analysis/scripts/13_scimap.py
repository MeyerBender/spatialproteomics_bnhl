import scimap as sm
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

# Read input from Snakemake
adata = sc.read_h5ad(snakemake.input.anndata_path)
adata.obs_names_make_unique()
label_col = snakemake.wildcards.label  # Get label from wildcard

# === Spatial distance ===
adata = sm.tl.spatial_distance(adata, x_coordinate='centroid-1', y_coordinate='centroid-0', phenotype=label_col, imageid='sample_id')

# === Spatial interaction ===
adata = sm.tl.spatial_interaction(
    adata, 
    x_coordinate='centroid-1', y_coordinate='centroid-0', phenotype=label_col,
    method='radius', 
    radius=140, 
    permutation=100,
    imageid='sample_id',
    label='spatial_interaction_radius'
)

# Write output
adata.write(snakemake.output.anndata_path)
