import xarray as xr
import spatialproteomics as sp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
import gc
import yaml

metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)
files = snakemake.input.zarr_files

# === reading the parameter config ===
with open(snakemake.input.parameter_config_path, "r") as yaml_file:
    yaml_data = yaml.safe_load(yaml_file)
k = yaml_data["k"]
r = yaml_data["radius"]

# === reading the zarrs into the sample dict ===
sample_dict = {}
for i, (sample_id_a, sample_id_b) in metadata[["sample_id_a", "sample_id_b"]].iterrows():
    # adding sample a and b to the dictionary if we have images for them
    for sample_id in [sample_id_a, sample_id_b]:
        if pd.isna(sample_id):
            continue
        file_mask = [f"{sample_id}_" in x for x in files]
        matching_files = [files[i] for i, mask in enumerate(file_mask) if mask]
        if len(matching_files) != 1:
            raise ValueError(f"Matching files for sample {sample_id}: {matching_files}")
        else:
            file = matching_files[0]
            sample_dict[sample_id] = xr.open_zarr(file)

image_container = sp.ImageContainer(sample_dict)
sample_dict = image_container.compute_neighborhoods(neighborhood_method='radius', radius=r, k=k, overwrite=False, seed=0)

# === computing heterogeneity both per cells and per samples ===
# Initialize an empty list to collect data for each sample
rows = []

# Iterate over samples and extract graph features
for sample in tqdm(sample_dict.keys()):
    ds = sample_dict[sample].nh.add_neighborhood_obs(
        features=['degree', 'homophily', 'inter_label_connectivity', 'diversity_index']
    )
    
    ds.drop_encoding().to_zarr(f"{snakemake.output.zarr_path}/{sample}.zarr")
    
    graph_features = ds.nh.compute_graph_features(
        features=['num_nodes', 'num_edges', 'density', 'modularity', 'assortativity']
    )
    
    # Add the sample name to the graph features dictionary
    graph_features['sample'] = sample
    
    # Append the dictionary to the rows list
    rows.append(graph_features)
    
    gc.collect()

# Convert the list of dictionaries to a DataFrame
graph_features_df = pd.DataFrame(rows)
graph_features_df.to_csv(snakemake.output.network_feature_df)
