import xarray as xr
import spatialproteomics as sp
import pandas as pd
import numpy as np
from scipy.ndimage import zoom

def upsample_mask(mask, upsample_factor, target_shape):
    # Upsample using nearest-neighbor interpolation
    upsampled = zoom(mask, zoom=upsample_factor, order=0)
    
    # Get current and target shapes
    up_h, up_w = upsampled.shape
    tgt_h, tgt_w = target_shape

    # Crop if larger
    upsampled = upsampled[:tgt_h, :tgt_w]

    # Pad if smaller
    pad_h = tgt_h - upsampled.shape[0]
    pad_w = tgt_w - upsampled.shape[1]

    if pad_h > 0 or pad_w > 0:
        upsampled = np.pad(
            upsampled,
            ((0, pad_h), (0, pad_w)),
            mode='constant',
            constant_values=0  # Adjust padding value if needed
        )

    return upsampled

zarr_path = snakemake.input.zarr_path
threshold_path = snakemake.input.threshold_path
mask_path = snakemake.input.mask_path

artefact_mask = np.load(mask_path)

# basic formatting of sample IDs etc.
file_name = zarr_path.split('/')[-1]
sample_id = file_name.split('/')[-1]
sample_id = f"{sample_id.split('_')[2]}_{sample_id.split('_')[3]}_{sample_id.split('_')[4]}".replace('_Scan1.zarr', '')

# reading in the files
ds = xr.open_zarr(zarr_path)
threshold_df = pd.read_csv(threshold_path, index_col=0)

upsampled_mask = upsample_mask(artefact_mask, upsample_factor=16, target_shape=(ds.dims['y'], ds.dims['x']))
assert upsampled_mask.shape == (ds.dims['y'], ds.dims['x']),f"Found shape {upsampled_mask.shape} for mask, but image had shape {(ds.dims['y'], ds.dims['x'])}."

# masking out the problematic cells and removing any outlying cells
ds = ds.pp.add_layer(upsampled_mask).pp.mask_cells().pp.remove_outlying_cells()

# copying the raw image and putting it into its own slot
img_raw = ds['_image'].values
ds = ds.pp.add_layer(img_raw, key_added='_image_raw')

# running thresholding
markers = threshold_df['marker'].values
thresholds = threshold_df['threshold'].values
ds = ds.pp.threshold(thresholds, channels=markers)

# performing quantification and running ct prediction
ct_marker_dict = {'PAX5': 'B cell', 'CD3': 'T cell', 'CD11b': 'Myeloid cell', 'CD11c': 'Dendritic cell', 'CD68': 'Macrophage', 'CD90': 'Stromal cell', 'Podoplanin': 'Stromal cell', 'CD31': 'Endothelial cell', 'CD34': 'Endothelial cell'}
ds = ds.pp.add_quantification(func='intensity_mean').pp.transform_expression_matrix(method='arcsinh').la.predict_cell_types_argmax(ct_marker_dict)

# write to zarr
ds.drop_encoding().to_zarr(snakemake.output.zarr_path)
