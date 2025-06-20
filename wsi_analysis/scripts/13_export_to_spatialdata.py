import xarray as xr
import spatialproteomics as sp
import spatialdata as sd
import tifffile

ds = xr.open_zarr(snakemake.input.zarr_path)
ds["_image"] = ds["_image"].astype("uint8")
markers = ds.coords['channels'].values
seg_stardist = tifffile.imread(snakemake.input.stardist_path)
sdata = ds.tl.convert_to_spatialdata()
sdata.images['image_raw'] = sd.models.Image2DModel.parse(ds['_image_raw'].values, transformations=None, dims=("c", "x", "y"), c_coords=markers)
print(ds['_image_raw'].values.shape)
print(seg_stardist.shape)
sdata.labels['segmentation_stardist'] = sd.models.Labels2DModel.parse(seg_stardist, transformations=None, dims=("x", "y"))
sdata.write(snakemake.output.spatialdata_path)