import xarray as xr
import spatialproteomics as sp
import pandas as pd
import numpy as np

zarr_path = snakemake.input.zarr_path
threshold_path = snakemake.input.threshold_path

# basic formatting of sample IDs etc.
sample_id = zarr_path.split('/')[-1].replace('.zarr', '')
print(sample_id)

# reading in the files
ds = xr.open_zarr(zarr_path)
threshold_df = pd.read_csv(threshold_path, index_col=0)

# runing thresholding
threshold_df = threshold_df[threshold_df['sample_id'] == sample_id]

markers_for_binarization = ['TCF7/TCF1' if marker == 'TCF7-TCF1' else marker for marker in threshold_df['marker'].values]
thresholds_for_binarization = threshold_df['quantile'].values
ds = ds.pp.threshold(intensity=thresholds_for_binarization, channels=markers_for_binarization)
obj_binarized = ds.pp.add_quantification(func=sp.percentage_positive, key_added='_perc_pos_filtered')

# going through every marker for which we have thresholds
binarization_results = dict()
for i, row in threshold_df.iterrows():
    # the original naming is confusion, threshold actually refers to the percentage of positive pixels, this should clear things up a bit
    marker, threshold, perc_pos =  row['marker'], row['quantile'], row['perc_pos']
    marker_name = 'TCF7/TCF1' if marker == 'TCF7-TCF1' else marker

    # Filter by lower value
    # in case the highest value is above the highest raw value, we set the binarization to 0
    try:
        obj_binarized = ds.pp[marker_name].pp.add_quantification(func=sp.percentage_positive, key_added='_perc_pos_filtered')
        # Binarize with the thresholds
        obj_binarized = obj_binarized.la.threshold_labels({marker_name: perc_pos}, layer_key='_perc_pos_filtered')
        # Extracting the binarization information
        binarized_cells = obj_binarized.pp.get_layer_as_df("_obs")[f"{marker_name}_binarized"]
    except AssertionError:
        binarized_cells = 0
    binarization_results[f"{marker}_binarized"] = binarized_cells

binarization_results = pd.DataFrame(binarization_results)

# === calling cell subtypes ===
# adding the binarization into the obs table
print(binarization_results.shape)
print(ds.sizes)
print(len(binarization_results))
print(binarization_results.head())
ds_final = ds.pp.add_obs_from_dataframe(binarization_results)

# predicting cell subtypes
ds_final = ds_final.la.predict_cell_subtypes(snakemake.input.gating_tree)

# write to zarr
ds_final.drop_encoding().to_zarr(snakemake.output.zarr_path)
