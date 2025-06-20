import pandas as pd

df_paths = snakemake.input.co_occurrence_dfs

df_full = []
for path in df_paths:
    df_full.append(pd.read_csv(path, index_col=0))
    
df_full = pd.concat(df_full, axis=0).reset_index(drop=True)
df_full.to_csv(snakemake.output.co_occurrence_df)