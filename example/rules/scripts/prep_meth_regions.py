import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BSSAMPLES = [
    'NPC_9sca_ctrl_BS_rep1',
    'NPC_9sca_ctrl_BS_rep2',
    'NPC_9sca_ko_BS_rep1',
    'NPC_9sca_ko_BS_rep2',
]
UPSTREAM = 3000
DOWNSTREAM = 3000
MDIFF = 10
LOW_METH = 10

def load_gene_set(p):
    path = Path(p)
    return set(pd.read_table(path, header=0, index_col=0).index)

groups = {
    'up': load_gene_set(snakemake.input.up),
    'down': load_gene_set(snakemake.input.down),
    'nonde': load_gene_set(snakemake.input.nonde),
}
gene_to_group = {g: grp for grp, genes in groups.items() for g in genes}

windows = []
with open(snakemake.input.gtf) as f:
    for line in f:
        if line.startswith('#'):
            continue
        li = line.rstrip('\n').split('\t')
        if li[2] != 'gene':
            continue
        gid = next((a.strip().split(' ')[1].strip('"') for a in li[8].split(';') if 'gene_id' in a), None)
        if gid not in gene_to_group:
            continue
        strand = li[6]
        tss = int(li[3]) if strand == '+' else int(li[4])
        wstart, wend = (tss - UPSTREAM, tss + DOWNSTREAM) if strand == '+' else (tss - DOWNSTREAM, tss + UPSTREAM)
        windows.append((li[0], max(wstart, 0), wend, gid, gene_to_group[gid]))

windows = pd.DataFrame(windows, columns=['chrom', 'start', 'end', 'gene_id', 'group'])
windows['start'] = windows['start'].astype('int32')
windows['end'] = windows['end'].astype('int32')
print(len(windows), 'gene windows')

def load_bedgraph(path, chroms):
    df = pd.read_csv(
        path, sep='\t', skiprows=1, header=None,
        names=['chrom', 'start', 'end', 'pct'],
        dtype={'chrom': str, 'start': 'int32', 'end': 'int32', 'pct': 'float32'},
    )
    return df[df['chrom'].isin(chroms)]

def sites_in_windows(sites, windows):
    results = []
    for chrom, w in windows.groupby('chrom', sort=False):
        s = sites[sites['chrom'] == chrom]
        if s.empty:
            continue
        w = w.sort_values('start')
        w_start = w['start'].to_numpy()
        w_end = w['end'].to_numpy()
        gene_id = w['gene_id'].to_numpy()
        group = w['group'].to_numpy()

        s_start = s['start'].to_numpy()
        idx = np.searchsorted(w_start, s_start, side='right') - 1
        valid = idx >= 0
        within = np.zeros(len(s), dtype=bool)
        within[valid] = s_start[valid] < w_end[idx[valid]]

        hit = s.loc[within].copy()
        hit_idx = idx[within]
        hit['gene_id'] = gene_id[hit_idx]
        hit['group'] = group[hit_idx]
        results.append(hit[['chrom', 'start', 'gene_id', 'group', 'pct']])
    if not results:
        return pd.DataFrame(columns=['chrom', 'start', 'gene_id', 'group', 'pct'])
    return pd.concat(results, ignore_index=True)

chroms = set(windows['chrom'])
per_sample = {}
for sample, bg_path in zip(BSSAMPLES, snakemake.input.bgs):
    bg = load_bedgraph(bg_path, chroms)
    hits = sites_in_windows(bg, windows)
    per_sample[sample] = hits.set_index(['chrom', 'start', 'gene_id', 'group'])['pct'].rename(sample)
    print(sample, len(hits), 'sites in gene windows')

sites = pd.concat(per_sample.values(), axis=1, join='outer').reset_index()
ctrl_cols = [c for c in BSSAMPLES if '_ctrl_' in c]
ko_cols = [c for c in BSSAMPLES if '_ko_' in c]
sites = sites.dropna(subset=ctrl_cols, how='all').dropna(subset=ko_cols, how='all')
sites['ctrl_mean'] = sites[ctrl_cols].mean(axis=1)
sites['ko_mean'] = sites[ko_cols].mean(axis=1)
sites['diff'] = sites['ko_mean'] - sites['ctrl_mean']

def top_site_demeth_ko(df):
    # 'up' = higher expression in ko (DE.R: logFC is ko relative to ctrl) ->
    # expect the ko sample specifically demethylated at this gene.
    idx = df.groupby('gene_id')['diff'].idxmin()
    out = df.loc[idx]
    return out[out['diff'] <= -MDIFF]

def top_site_demeth_ctrl(df):
    # 'down' = higher expression in ctrl -> expect the ctrl sample
    # specifically demethylated at this gene.
    idx = df.groupby('gene_id')['diff'].idxmax()
    out = df.loc[idx]
    return out[out['diff'] >= MDIFF]

def top_site_unmeth(df):
    low = df[(df['ctrl_mean'] <= LOW_METH) & (df['ko_mean'] <= LOW_METH)]
    score = low[['ctrl_mean', 'ko_mean']].max(axis=1)
    idx = score.groupby(low['gene_id']).idxmin()
    return low.loc[idx]

selected = pd.concat([
    top_site_demeth_ko(sites[sites['group'] == 'up']),
    top_site_demeth_ctrl(sites[sites['group'] == 'down']),
    top_site_unmeth(sites[sites['group'] == 'nonde']),
], ignore_index=True)

print(selected.groupby('group').size())
selected.groupby('group')[['ctrl_mean', 'ko_mean', 'diff']].describe()

for of, grp in zip([snakemake.output.up, snakemake.output.down, snakemake.output.nonde], ['up', 'down', 'nonde']):
    sub = selected[selected['group'] == grp].sort_values(['chrom', 'start'])
    bed = sub[['chrom', 'start', 'gene_id']].copy()
    bed['end'] = sub['start'] + 1
    bed = bed[['chrom', 'start', 'end']]
    bed.to_csv(of, sep='\t', header=False, index=False)
