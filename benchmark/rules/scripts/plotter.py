import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

_stem_re = re.compile(r'_(dt[34])_t(\d+)$')

def read_benchmark(file):
    stem = file.split('/')[-1].split('.')[0]
    m = _stem_re.search(stem)
    a = pd.read_table(file, sep='\t')
    a = a[['s', 'max_rss']]
    a['rep'] = a.index + 1
    a['version'] = m.group(1)
    a['threads'] = int(m.group(2))
    a['sample'] = stem[:m.start()]
    a['os'] = snakemake.params.os
    a['modality'] = file.split('/')[-2]
    return a

_dfs = []

for file in snakemake.input:
    _dfs.append(read_benchmark(file))

a = pd.concat(_dfs, ignore_index=True)
a.to_csv(snakemake.output.csv, sep='\t', index=False)


# Plot figure (max threads per modality, i.e. non-exhaustive comparison)
b = a[a['threads'] == a.groupby('modality')['threads'].transform('max')]

_mods = b['modality'].unique()
ap = b.pivot(index=['rep', 'sample', 'os', 'modality'], columns='version', values=['s', 'max_rss']).reset_index()
ap['speedup'] = ap['s']['dt3']/ap['s']['dt4']
ap['memory'] = ap['max_rss']['dt3']/ap['max_rss']['dt4']

fig, ax = plt.subplots(nrows=len(_mods), ncols=3, figsize=(24,16))

_axix = 0
for mod in _mods:
    dt3 = b[(b['modality'] == mod) & (b['version'] == 'dt3')]
    dt4 = b[(b['modality'] == mod) & (b['version'] == 'dt4')]

    ax[_axix, 0].plot(
        ['v3','v4'],
        [dt3['s'], dt4['s']]
    )
    ax[_axix, 1].plot(
        ['v3','v4'],
        [dt3['max_rss'], dt4['max_rss']]
    )

    sns.scatterplot(
        data=ap[ap['modality'] == mod],
        x='speedup',
        y='memory',
        ax=ax[_axix, 2]
    )
    ax[_axix, 2].set_xlabel('speedup (time)')
    ax[_axix, 2].set_ylabel('memory dt3 / memory dt4')

    ax[_axix, 0].set_title(f"{mod} - time")
    ax[_axix, 1].set_title(f"{mod} - memory")
    _axix += 1

fig.savefig(snakemake.output.png, dpi=300)
