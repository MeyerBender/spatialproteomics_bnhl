import tifffile 
import spatialproteomics as sp
import pandas as pd
import yaml

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path
segmentation_path = snakemake.input.mask_path
parameter_config_path = snakemake.input.parameter_config_path

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers)

# === parameter config ===
with open(parameter_config_path, "r") as yaml_file:
    yaml_data = yaml.safe_load(yaml_file)
mask_growth = yaml_data["segmentation"]["mask_growth_pixels"]
min_cell_size = yaml_data["segmentation"]["min_cell_size"]
max_cell_size = yaml_data["segmentation"]["max_cell_size"]

# === fixing disconnected cells ===
segmentation = tifffile.imread(segmentation_path)
img = img.pp.add_segmentation(segmentation)

# === filtering segmentation ===
img = img.pp.add_observations("area").pp.filter_by_obs('area', func=lambda x: (x >= min_cell_size) & (x <= max_cell_size))

# === growing segmentation ===
img = img.pp.grow_cells(iterations=mask_growth)

# === exporting ===
segmentation_processed = img['_segmentation'].values
tifffile.imwrite(snakemake.output.mask_path, segmentation_processed)