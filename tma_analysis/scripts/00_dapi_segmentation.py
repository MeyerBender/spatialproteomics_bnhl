import tifffile 
import spatialproteomics as sp
import pandas as pd

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers)

# === segmentation ===
img = img.tl.cellpose(channels="DAPI")

# === filtering segmentation ===
img = img.pp.add_observations("area").pp.filter_by_obs('area', func=lambda x: (x >= 50) & (x <= 300))

# === growing segmentation ===
img = img.pp.grow_cells(iterations=mask_growth)
    
# === exporting ===
img.drop_encoding().to_zarr(output_path)
