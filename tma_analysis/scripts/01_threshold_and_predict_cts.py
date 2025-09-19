import xarray as xr
import spatialproteomics as sp
import numpy as np
import pandas as pd
import os
import tifffile

threshold_path = snakemake.input.threshold_path
zarr_path = snakemake.input.zarr_path
gating_tree_path = snakemake.input.gating_tree_path
output_path = snakemake.output.zarr_path

# this is not in snakemake, since not every sample has a mask
mask_path = '/g/huber/users/meyerben/data/codex/250128_bnhl_artefact_masks'

sample_id = zarr_path.split('/')[-1].replace('.zarr', '')
sample_prefix = '_'.join(sample_id.split('_')[:-1])

# === dictionary mapping from cell types to markers ===CD
# loading the thresholds from multirazor
marker_ct_dict = {'PAX5': 'B cell', 'CD3': 'T cell', 'CD11b': 'Myeloid cell', 'CD11c': 'Dendritic cell', 'CD68': 'Macrophage', 'CD90': 'Stromal cell', 'Podoplanin': 'Stromal cell', 'CD31': 'Endothelial cell', 'CD34': 'Endothelial cell'}
channels = list(marker_ct_dict.keys())
labels = list(marker_ct_dict.values())

# === loading the binarazor df ===
threshold_df = pd.read_csv(threshold_path)
# removing unnecessary channels (not doing this bc we might need the others for the binarization)
# threshold_df = threshold_df[threshold_df['channel'].isin(ct_channel_dict.values())]
# setting a prefix column to be in line with the metadata
threshold_df['prefix'] = ['_'.join(x.split('_')[:-1]) for x in threshold_df['sample']]
# only keeping samples which we have in the metadata
# threshold_df = threshold_df[threshold_df['prefix'].isin(metadata['prefix'].values)]
# only keeping thresholds which were reviewed
threshold_df = threshold_df[threshold_df['status'] == 'reviewed']

# === calling cell types ===
quantiles = []
valid_sample = True
for channel in channels:
    if len(threshold_df.loc[(threshold_df['prefix'] == sample_prefix) & (threshold_df['channel'] == channel), 'threshold']) == 1:
        quantile = threshold_df.loc[(threshold_df['prefix'] == sample_prefix) & (threshold_df['channel'] == channel), 'threshold'].iloc[0]
        quantiles.append(quantile)
    else:
        print(f"Could not find threshold for sample {sample_prefix}, channel {channel}, skipping sample!")
        valid_sample = False

# loading the dataset and removing the outlying cells
ds = xr.load_dataset(zarr_path, engine='zarr').pp.remove_outlying_cells()
        
if valid_sample:
    ds = (ds.pp.threshold(quantiles, channels=channels)
                .pp.add_quantification(func='intensity_mean')
                .pp.transform_expression_matrix(method='arcsinh')
                .la.predict_cell_types_argmax(marker_ct_dict)
                )

    # === binarizing functional markers ===
    # dropping NAs here because the exported table contains some channels for which we do not have valid thresholds
    threshold_df_tmp = threshold_df[(threshold_df['sample'] == sample_id) & (threshold_df['status'] == 'reviewed')].dropna()
    # Exclude DAPI
    threshold_df_tmp = threshold_df_tmp[threshold_df_tmp['channel'] != 'DAPI']

    markers_for_binarization = ['TCF7/TCF1' if marker == 'TCF7-TCF1' else marker for marker in threshold_df_tmp['channel'].values]
    thresholds_for_binarization = threshold_df_tmp['lower'].values
    ds = ds.pp.threshold(intensity=thresholds_for_binarization, channels=markers_for_binarization)
    obj_binarized = ds.pp.add_quantification(func=sp.percentage_positive, key_added='_perc_pos_filtered')

    # going through every marker for which we have thresholds
    binarization_results = dict()
    for i, row in threshold_df_tmp.iterrows():
        # the original naming is confusion, threshold actually refers to the percentage of positive pixels, this should clear things up a bit
        marker, threshold, perc_pos =  row['channel'], row['lower'], row['threshold']
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
    ds_final = ds.pp.add_obs_from_dataframe(binarization_results)
    
    # predicting cell subtypes
    ds_final = ds_final.la.predict_cell_subtypes(gating_tree_path)

    # === removing masked regions ===
    if os.path.exists(f"{mask_path}/{sample_id}.tiff"):
        print(f"Found mask for {sample_id}, adding that to the object")
        mask = tifffile.imread(f"{mask_path}/{sample_id}.tiff")
        # add the mask to the spatialproteomics object
        ds_final = ds_final.pp.add_layer(mask).pp.mask_cells()
else:
    # we still want to write out the unprocessed zarrs for files for which we do not have thresholds
    ds_final = ds
    
# === storing the resulting dataset ===
# redoing quantification so that the functional markers are also quantified properly
ds_final.pp.drop_layers('_intensity').pp.add_quantification(func='intensity_mean').pp.transform_expression_matrix(method='arcsinh').drop_encoding().to_zarr(output_path)