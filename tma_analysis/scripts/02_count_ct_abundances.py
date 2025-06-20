import pandas as pd
import numpy as np
import xarray as xr
import spatialproteomics as sp

metadata = pd.read_csv(snakemake.input.metadata_path, index_col=0)
files = snakemake.input.zarr_files

# === reading the zarrs into the sample dict ===
sample_dict = {}
for i, (sample_id_a, sample_id_b) in metadata[["sample_id_a", "sample_id_b"]].iterrows():
    # adding sample a and b to the dictionary if we have images for them
    for sample_id in [sample_id_a, sample_id_b]:
        if pd.isna(sample_id):
            continue
        file_mask = [f"{sample_id}_" in x for x in files]
        matching_files = [files[i] for i, mask in enumerate(file_mask) if mask]
        if len(matching_files) != 1:
            raise ValueError(f"Matching files for sample {sample_id}: {matching_files}")
        else:
            file = matching_files[0]
            sample_dict[sample_id] = xr.open_zarr(file)

# === computing absolute abundances ===
def compute_abundance_df(level, cts):
    abundance_df_per_patient = pd.DataFrame(np.zeros((metadata.shape[0], len(cts))), columns=cts)
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
                ct_predictions = sample_dict[codex_id].pp.get_layer_as_df('_la_layers')[f'labels_{level}'].value_counts()
                for ct in list(cts):
                    abundance_df_per_patient.loc[patient_id, ct] += ct_predictions.get(ct, 0)
            except AssertionError:
                print(f"Could not find thresholds for sample {codex_id}, continuing...")

    # remove rows with 0 of all cells (these are the unmapped patients)
    abundance_df_per_patient = abundance_df_per_patient[abundance_df_per_patient.sum(axis=1) != 0]
        
    print(abundance_df_per_patient.shape)
    return abundance_df_per_patient

# === storing the results ===
df = compute_abundance_df(level=0, cts=['B cell', 'T cell', 'Myeloid cell', 'Dendritic cell', 'Macrophage', 'Stromal cell', 'Endothelial cell'])
df.to_csv(snakemake.output.level_0_abundance_path)

df = compute_abundance_df(level=1, cts=['B cell', 'B_prol', 'T_h', 'T_tox', 'Myeloid cell', 'Dendritic cell', 'Macrophage', 'M2', 'Stromal cell', 'Endothelial cell'])
df.to_csv(snakemake.output.level_1_abundance_path)

df = compute_abundance_df(level=2, cts=['B cell', 'B_prol', 'T_h', 'T_h_mem', 'T_h_naive', 'T_reg', 'T_fh', 'T_tox', 'T_exh', 'T_tox_naive', 'T_tox_mem', 'Myeloid cell', 'Dendritic cell', 'Macrophage', 'M2', 'Stromal cell', 'Endothelial cell'])
df.to_csv(snakemake.output.level_2_abundance_path)

df = compute_abundance_df(level=3, cts=['B cell', 'B_prol', 'T_h', 'T_h_mem', 'T_h_naive', 'T_reg', 'T_fh', 'T_reg_Helios', 'T_tox', 'T_exh', 'T_progenitor_exh', 'T_terminally_exh', 'T_tox_naive', 'T_tox_mem', 'Myeloid cell', 'Dendritic cell', 'Macrophage', 'M2', 'Stromal cell', 'Endothelial cell'])
df.to_csv(snakemake.output.level_3_abundance_path)
