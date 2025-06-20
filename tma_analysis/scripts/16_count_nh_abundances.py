import pandas as pd
import numpy as np
import xarray as xr
import spatialproteomics as sp

metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)
files = snakemake.input.zarr_files

nh_mapping_dict = {'Neighborhood 0': 'T/Dendritic',
                      'Neighborhood 1': 'T',
                      'Neighborhood 2': 'Mixed',
                      'Neighborhood 3': 'B_prol/T_fh',
                      'Neighborhood 4': 'T_h/T_reg',
                      'Neighborhood 5': 'B',
                      'Neighborhood 6': 'Myeloid/Macro'}

neighborhood_colors = {'B': '#5799d1',
                       'B_prol/T_fh': '#154e91',
                       'T/Dendritic': '#fef1c3',
                       'T': '#ebc850',
                       'T_h/T_reg': '#cca02d',
                       'Myeloid/Macro': '#de6866',
                       'Mixed': '#c8ceda'}

# === reading the zarrs into the sample dict ===
sample_dict = {}
for i, (sample_id_a, sample_id_b) in metadata[["sample_id_a", "sample_id_b"]].iterrows():
    # adding sample a and b to the dictionary if we have images for them
    for sample_id in [sample_id_a, sample_id_b]:
        if pd.isna(sample_id):
            continue
        file_mask = [f"{sample_id}." in x for x in files]
        matching_files = [files[i] for i, mask in enumerate(file_mask) if mask]
        if len(matching_files) != 1:
            raise ValueError(f"Matching files for sample {sample_id}: {matching_files}")
        else:
            file = matching_files[0]
            sample_dict[sample_id] = xr.open_zarr(file)

# === computing absolute abundances ===
def compute_abundance_df(nhs):
    abundance_df_per_patient = pd.DataFrame(np.zeros((metadata.shape[0], len(nhs))), columns=nhs)
    abundance_df_per_patient.index = list(metadata['Histo-Nr'].values)

    for i, row in metadata.iterrows():
        patient_id = row['Histo-Nr']

        # getting the matching codex ids. One is already in prefix, the other requires a bit more manipulation
        matching_codex_ids = [row['sample_id_a'], row['sample_id_b']]

        if len(matching_codex_ids) == 0:
            print(f"Could not find matching CODEX image for {patient_id}, continuing...")
            continue

        # going through the matching codex ids and computing the required values
        abundance_dict = {}    
        for codex_id in matching_codex_ids:
            if codex_id not in sample_dict.keys():
                print(f"Could not find matching CODEX image for {patient_id}/{codex_id}, continuing...")
                continue
            try:
                nh_predictions = sample_dict[codex_id].pp.get_layer_as_df()['_neighborhoods']
                nh_predictions = pd.Series([nh_mapping_dict[x] for x in nh_predictions])
                ct_predictions = nh_predictions.value_counts()
                for nh in list(nhs):
                    abundance_df_per_patient.loc[patient_id, nh] += ct_predictions.get(nh, 0)
            except AssertionError:
                print(f"Could not find thresholds for sample {codex_id}, continuing...")

    # remove rows with 0 of all cells (these are the unmapped patients)
    abundance_df_per_patient = abundance_df_per_patient[abundance_df_per_patient.sum(axis=1) != 0]
        
    print(abundance_df_per_patient.shape)
    return abundance_df_per_patient

# === storing the results ===
df = compute_abundance_df(nhs=['B', 'B_prol/T_fh', 'Mixed', 'Myeloid/Macro', 'T', 'T/Dendritic', 'T_h/T_reg'])
df.to_csv(snakemake.output.abundance_path)
