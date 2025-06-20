import tifffile 
import spatialproteomics as sp
import pandas as pd
import gc

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path

diameters = {'20240923_DLBCL_1XFLMK_B_Scan1': 13.68, '20240920_DLBCL_E13CKP_5_Scan1': 13.63, '20240917_DLBCL_7SS1X1_IA3_2_Scan1': 16.66, '20240920_DLBCL_5MXJ3L_A_Scan1': 15.05, '20240923_DLBCL_Z979KU_3_Scan1': 17.00}

sample = image_path.split("/")[-1].replace(".qptiff", "")
diameter = diameters.get(sample, None)

print(f"Got diameter: {diameter}")

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers).pp["DAPI"]
gc.collect()

# === segmentation ===
img = img.tl.cellpose(channel="DAPI", diameter=diameter)
    
# === exporting ===
segmentation_raw = img['_cellpose_segmentation'].values
tifffile.imwrite(snakemake.output.mask_path, segmentation_raw)