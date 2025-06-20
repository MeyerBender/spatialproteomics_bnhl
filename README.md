# Analysing B cell Non-Hodgkin lymphomas with spatialproteomics
This repository contains the code used to process the data and create the figures for the manuscript "Spatialproteomics: an interoperable toolbox for analyzing highly multiplexed fluorescence image data" ([manuscript](https://www.biorxiv.org/content/10.1101/2025.04.29.651202v1)).
To see how you can apply `spatialproteomics` to your own data, please refer to the [official documentation](https://haraldvohringer.com/spatialproteomics/notebooks/ExampleWorkflow.html).

### TMA analysis
The folder `tma_analysis` contains the snakemake pipeline which was used to process the original tissue microarray images, from qptiff to zarr. The data derived from this pipeline can be found on the BioImage archive (link coming soon). The structure of this pipeline is highly dependent on the available computing infrastructure. For a more detailed tutorial on how to use `spatialproteomics` in combination with `snakemake` for your own data, please refer to [this repository](https://github.com/MeyerBender/spatialproteomics_pipeline_template).

### WSI analysis
The folder `wsi_analysis` contains the snakemake pipeline which was used to process the original whole slide images, from qptiff to zarr. The data derived from this pipeline can be found on the BioImage archive.

### Figures
This folder contains the code used to generate the figures in the manuscript.
