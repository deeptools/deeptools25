import pandas as pd
import seaborn as sns
from pathlib import Path

def read_benchmark(file):
    _version = file.split('/')[-1].split('.')[0].split('_')[-1]
    sample = file.split('/')[-1].split('.')[0].replace('_dt3', '').replace('_dt4', '')
    a = pd.read_table(file, sep='\t')    
    a = a[['s', 'max_rss']]
    a['rep'] = a.index + 1
    a['version'] = _version
    a['sample'] = sample
    a['os'] = snakemake.params.os
    a['modality'] = file.split('/')[-2]
    return a

_dfs = []

for file in snakemake.input:
    _dfs.append(read_benchmark(file))

df = pd.concat(_dfs, ignore_index=True)
print(df)

df.to_csv(snakemake.output.csv, sep='\t', index=False)
Path(snakemake.output.png).touch()