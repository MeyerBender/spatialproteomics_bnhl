import xarray as xr
import spatialproteomics as sp
import pandas as pd
import numpy as np

ds = xr.open_zarr(snakemake.input.zarr_path)
sample_id = snakemake.input.zarr_path.split('/')[-1].replace('.zarr', '')
sample_id_list = sample_id.split('_')
sample_id = f'{sample_id_list[0]}_{sample_id_list[1]}_{sample_id_list[2]}'

df1 = ds.pp.get_layer_as_df()[['centroid-0', 'centroid-1']]
try:
    df2 = ds.pp.get_layer_as_df('_la_layers')
    df = pd.concat((df1, df2), axis=1)
    df['sample'] = sample_id
    df.to_csv(snakemake.output.csv_path, index=False)
except AssertionError:
    # if we have not called any cell types because of QC issues
    df1.to_csv(snakemake.output.csv_path, index=False)
