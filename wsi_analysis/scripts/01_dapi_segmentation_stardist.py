import tifffile 
import spatialproteomics as sp
import pandas as pd
import gc

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path

sample = image_path.split("/")[-1].replace(".qptiff", "")

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers).pp["DAPI"]
gc.collect()

# === segmentation ===
img = img.tl.stardist(channel="DAPI")
    
# === exporting ===
segmentation_raw = img['_stardist_segmentation'].values
tifffile.imwrite(snakemake.output.mask_path, segmentation_raw)