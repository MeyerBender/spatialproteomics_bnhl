import tifffile 
import spatialproteomics as sp
import pandas as pd
import yaml

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path
segmentation_path = snakemake.input.segmentation_path

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers)
segmentation = tifffile.imread(segmentation_path)
img = img.pp.add_segmentation(segmentation)

# === exporting ===
img.drop_encoding().to_zarr(snakemake.output.zarr_path)