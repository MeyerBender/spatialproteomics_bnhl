import xarray as xr
import spatialproteomics as sp
import anndata as ad

ds = xr.open_zarr(snakemake.input.zarr_path)
adata = ds.tl.convert_to_anndata()

# also getting the lower level ct predictions
df = ds.pp.get_layer_as_df('_la_layers')
for label_col in df.columns:
    adata.obs[label_col] = df[label_col].values

adata.write(snakemake.output.adata_path)