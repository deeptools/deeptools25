import pandas as pd
import seaborn as sns

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

a = pd.concat(_dfs, ignore_index=True)
print(df)

a.to_csv(snakemake.output.csv, sep='\t', index=False)


# Plot figure
_mods = a['modality'].unique()
ap = a.pivot(index=['rep', 'sample', 'os', 'modality'], columns='version', values=['s', 'max_rss']).reset_index()
ap['speedup'] = ap['s']['dt3']/ap['s']['dt4']
ap['memory'] = ap['max_rss']['dt3']/ap['max_rss']['dt4']

fig, ax = plt.subplots(nrows=len(_mods), ncols=3, figsize=(24,16))

_axix = 0
for mod in _mods:
    dt3 = a[(a['modality'] == mod) & (a['version'] == 'dt3')]
    dt4 = a[(a['modality'] == mod) & (a['version'] == 'dt4')]
                
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
