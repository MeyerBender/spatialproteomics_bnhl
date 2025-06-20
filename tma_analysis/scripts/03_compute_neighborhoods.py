import xarray as xr
import spatialproteomics as sp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)
files = snakemake.input.zarr_files
k = int(snakemake.params.k)
r = int(snakemake.params.r)

# === reading the zarrs into the sample dict ===
sample_dict = {}
for i, (sample_id_a, sample_id_b) in tqdm(metadata[["sample_id_a", "sample_id_b"]].iterrows()):
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

samples = []
cells = []

for sample in image_container.objects.keys():
    cells_tmp = image_container.objects[sample].coords['cells'].values
    num_cells = len(cells_tmp)

    samples.extend([sample] * num_cells)
    cells.extend(list(cells_tmp))

cell_df = pd.DataFrame({'cell': cells, 'sample': samples})

df = pd.concat([cell_df, image_container.neighborhood_df.reset_index(drop=True), image_container.kmeans_df.reset_index(drop=True)], axis=1)

df.to_csv(snakemake.output.neighborhoods_csv)

# === plotting ===
# obtaining a data frame containing the neighborhood compositions
nh_composition = image_container.get_neighborhood_composition()
nh_composition.to_csv(snakemake.output.nh_composition_csv)

# plotting the neighborhood composition as a seaborn clustermap
cluster_map = sns.clustermap(
    nh_composition,
    cmap='coolwarm',
    cbar_pos=(1.02, 0.2, 0.03, 0.77),  # Position the colorbar to the right of the heatmap
    dendrogram_ratio=(0.00001, 0),  # Remove row dendrogram
    figsize=(9, 7),  # Adjust the size
    annot=False,  # No annotation
    center=0  # Setting the center of the colorbar to 0
)

# customizations of the plot
cluster_map.ax_heatmap.set_ylabel('')
plt.savefig(snakemake.output.plot_path)
