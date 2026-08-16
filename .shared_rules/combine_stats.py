import pandas as pd

df = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input], ignore_index=True)
df.to_csv(snakemake.output.tsv, sep="\t", index=False)
