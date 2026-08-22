import pandas as pd
import matplotlib.pyplot as plt


def parse_subreg(p):
    ixs = []
    with open(p) as f:
        for line in f:
            ix = line.strip().split('\t')[0]
            if ix != 'logFC':
                ixs.append(ix)
    return ixs


a = pd.read_table(snakemake.input.res, sep='\t', index_col=0)

up = parse_subreg(snakemake.input.up)
down = parse_subreg(snakemake.input.down)
nonde = parse_subreg(snakemake.input.nonde)

a['cat'] = 'ns'
a.loc[up, 'cat'] = 'up'
a.loc[down, 'cat'] = 'down'
a.loc[nonde, 'cat'] = 'ns - viz'
a['cat'] = pd.Categorical(a['cat'], categories=['ns', 'ns - viz', 'down', 'up'])

colors = {'ns': '#e0e0e0', 'ns - viz': '#8c8c8c', 'down': '#2166ac', 'up': '#b2182b'}
zorders = {'ns': 1, 'ns - viz': 2, 'down': 3, 'up': 3}

fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
for cat in a['cat'].cat.categories:
    sub = a[a['cat'] == cat]
    ax.scatter(sub['logCPM'], sub['logFC'], c=colors[cat], s=10, alpha=0.7,
               zorder=zorders[cat], label=cat)

ax.axhline(0, color='#e0e0e0', lw=0.8, linestyle='--', zorder=0)
ax.set_xlabel('logCPM')
ax.set_ylabel('logFC')

n_up, n_down, n_nsviz = len(up), len(down), len(nonde)

ax.text(0.98, 0.975, f'up: {n_up}', transform=ax.transAxes,
        ha='right', va='top', fontsize=10, color=colors['up'])
ax.text(0.98, 0.025, f'down: {n_down}', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=10, color=colors['down'])
ax.text(0.98, 0.55, f'ns - viz: {n_nsviz}', transform=ax.transAxes,
        ha='right', va='center', fontsize=10, color=colors['ns - viz'])
ax.set_xlim(-6, 20)

fig.savefig(snakemake.output.png, dpi=300)
