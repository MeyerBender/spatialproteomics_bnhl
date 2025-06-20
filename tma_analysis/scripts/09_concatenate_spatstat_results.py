import pandas as pd

pd.concat([pd.read_csv(f) for f in snakemake.input.ripley_csv]).to_csv(snakemake.output.ripley_merge)
pd.concat([pd.read_csv(f) for f in snakemake.input.pcf_csv]).to_csv(snakemake.output.pcf_merge)
pd.concat([pd.read_csv(f) for f in snakemake.input.cross_csv]).to_csv(snakemake.output.cross_merge)
pd.concat([pd.read_csv(f) for f in snakemake.input.pcf_cross_csv]).to_csv(snakemake.output.pcf_cross_merge)