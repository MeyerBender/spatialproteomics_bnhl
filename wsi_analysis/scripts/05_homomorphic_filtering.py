from glob import glob
import xarray as xr
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cv2
import spatialproteomics as sp

def homomorphic(img, radius = 13, blur=47):
    hh, ww = img.shape[:2]
    # take ln of image
    img_log = np.log1p(np.float64(img), dtype=np.float64)

    # do dft saving as complex output
    dft = np.fft.fft2(img_log, axes=(0,1))

    # apply shift of origin to center of image
    dft_shift = np.fft.fftshift(dft)

    # create black circle on white background for high pass filter
    mask = np.zeros_like(img, dtype=np.float64)
    cy = mask.shape[0] // 2
    cx = mask.shape[1] // 2
    cv2.circle(mask, (cx,cy), radius, 1, -1)
    mask = 1 - mask

    # antialias mask via blurring
    mask = cv2.GaussianBlur(mask, (blur, blur), 0)

    # apply mask to dft_shift
    dft_shift_filtered = np.multiply(dft_shift,mask)

    # shift origin from center to upper left corner
    back_ishift = np.fft.ifftshift(dft_shift_filtered)

    # do idft saving as complex
    img_back = np.fft.ifft2(back_ishift, axes=(0,1))

    # combine complex real and imaginary components to form (the magnitude for) the original image again
    img_back = np.abs(img_back)

    # apply exp to reverse the earlier log
    img_homomorphic = np.exp(img_back, dtype=np.float64)

    # scale result
    img_homomorphic = cv2.normalize(img_homomorphic, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    
    return img_homomorphic

zarr_path = snakemake.input.zarr_path

# basic formatting of sample IDs etc.
file_name = zarr_path.split('/')[-1]
sample_id = file_name.split('/')[-1]
sample_id = f"{sample_id.split('_')[2]}_{sample_id.split('_')[3]}_{sample_id.split('_')[4]}".replace('_Scan1.zarr', '')

# reading in the files
ds = xr.open_zarr(zarr_path)

# copying the raw image and putting it into its own slot
img_raw = ds['_image'].values
ds = ds.pp.add_layer(img_raw, key_added='_image_raw')

ds = ds.pp.apply(homomorphic)

# write to zarr
ds.drop_encoding().to_zarr(snakemake.output.zarr_path)
