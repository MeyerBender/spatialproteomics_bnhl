import anndata as ad
import scanpy as sc
import os
import pandas as pd
from tqdm.auto import tqdm
from glob import glob

# Load metadata
metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)

adatas = []  # List to store individual AnnData objects

for path in tqdm(snakemake.input.anndata_paths):
    # Read the AnnData file
    adata = sc.read_h5ad(path)
    
    # Extract sample ID from filename (assuming filename is sample-specific)
    sample_id = os.path.basename(path).replace(".h5ad", "")  # Remove extension
    sample_id = sample_id.rsplit('_', 1)[0]
    print(sample_id)
    
    # Find matching metadata row
    matching_rows = metadata[(metadata["sample_id_a"] == sample_id) | (metadata["sample_id_b"] == sample_id)]
    
    # Ensure exactly one match is found
    if len(matching_rows) != 1:
        raise ValueError(f"Expected exactly one metadata row for sample_id '{sample_id}', but found {len(matching_rows)}")

    # Get metadata as a dictionary (excluding sample_id_a and sample_id_b)
    metadata_dict = matching_rows.drop(columns=["sample_id_a", "sample_id_b"]).iloc[0].to_dict()
    
    # Store sample ID and metadata in obs
    adata.obs["sample_id"] = sample_id
    for key, value in metadata_dict.items():
        adata.obs[key] = value

    # Append to list
    adatas.append(adata)

# Concatenate all AnnData objects along observations (cells)
adata_combined = ad.concat(adatas, join="outer", label="sample_id", keys=[a.obs["sample_id"][0] for a in adatas])
adata_combined.write(snakemake.output.anndata_path)
