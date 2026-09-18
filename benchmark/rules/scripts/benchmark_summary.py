import re
import pandas as pd

_stem_re = re.compile(r'_(dt[34])_t(\d+)_rep(\d+)$')

def read_benchmark(file):
    stem = file.split('/')[-1].split('.')[0]
    m = _stem_re.search(stem)
    a = pd.read_table(file, sep='\t')
    a = a[['s', 'max_rss']]
    a['rep'] = int(m.group(3))
    a['version'] = m.group(1)
    a['threads'] = int(m.group(2))
    a['sample'] = stem[:m.start()]
    a['tool'] = file.split('/')[-2]
    return a

df = pd.concat([read_benchmark(f) for f in snakemake.input], ignore_index=True)
max_threads = snakemake.params.max_threads

# memory usage change (dt4 vs dt3), averaged over all threads/reps per sample
mem = (
    df.groupby(['tool', 'sample', 'version'])['max_rss']
      .mean()
      .unstack('version')
      .dropna(subset=['dt3', 'dt4'])
)
mem['mem_diff_gb'] = (mem['dt4'] - mem['dt3']) / 1024
mem['mem_diff_pct'] = 100 * (mem['dt4'] - mem['dt3']) / mem['dt3']
mem_summary = mem.groupby('tool')[['mem_diff_pct', 'mem_diff_gb']].mean()
mem_summary.columns = ['mean_mem_diff_pct', 'mean_mem_diff_gb']

# runtime speedup (dt3 / dt4) at max_threads, averaged per sample
t_max = df[df['threads'] == max_threads]
time = (
    t_max.groupby(['tool', 'sample', 'version'])['s']
         .mean()
         .unstack('version')
         .dropna(subset=['dt3', 'dt4'])
)
time['speedup'] = time['dt3'] / time['dt4']
speed_summary = time.groupby('tool')['speedup'].agg(['mean', 'median'])
speed_summary.columns = ['mean_speedup', 'median_speedup']

summary = mem_summary.join(speed_summary, how='outer').reset_index()
summary = summary[['tool', 'mean_speedup', 'median_speedup', 'mean_mem_diff_pct', 'mean_mem_diff_gb']]
summary = summary.round(2)
summary.to_csv(snakemake.output.tsv, sep='\t', index=False)
