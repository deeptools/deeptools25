import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
    a['os'] = snakemake.params.os
    a['modality'] = file.split('/')[-2]
    return a

_dfs = []

for file in snakemake.input:
    _dfs.append(read_benchmark(file))

df = pd.concat(_dfs, ignore_index=True)
df.to_csv(snakemake.output.csv, sep='\t', index=False)

df['organism'] = [i.split("_")[1] for i in df['sample']]
df['datatype'] = df['sample'].str.rsplit('_', n=1).str[1]
df['version'] = df['version'].map({'dt3': 'v3', 'dt4': 'v4'})
modalities = df["modality"].unique()
organisms  = ["human", "triticum"]
baseline = (
    df[df.version == 'v3']
    .groupby(['sample', 'threads'])[['s', 'max_rss']]
    .mean()
    .rename(columns={'s': 's_base', 'max_rss': 'rss_base'})
)
df = df.merge(baseline, on=['sample', 'threads'], how='left')
df['s_speedup']   = df['s_base'] / df['s']
df['rss_speedup'] = df['rss_base'] / df['max_rss']
def sem(x):
    return np.std(x, ddof=1) / np.sqrt(len(x))
def summarize(sub, cols):
    return (
        sub.groupby("threads")[cols]
           .agg(['mean', sem])
           .reset_index()
    )
org_colors    = {"human": "#5B8FC7", "triticum": "#C76B95"}
version_style = {"v3": dict(linestyle="--", marker="s"),
                  "v4": dict(linestyle="-",  marker="o")}
n_panels = len(modalities) + 1
ncols = -(-n_panels // 2)
fig, axes = plt.subplots(
    nrows=2,
    ncols=ncols,
    figsize=(12, 6),
    squeeze=False,
    constrained_layout=True
)
axes = axes.flatten()
for i, mod in enumerate(modalities):
    ax = axes[i]
    sub_mod = df[df["modality"] == mod]
    label_points = []
    for org in organisms:
        sub_org = sub_mod[sub_mod["organism"] == org]
        color = org_colors[org]
        last_x, last_y = None, None
        for v in ["v3", "v4"]:
            sub = sub_org[sub_org["version"] == v]
            if sub.empty:
                continue
            agg = summarize(sub, ["s_speedup", "rss_speedup"])
            threads = agg["threads"]
            ax.errorbar(
                threads, agg[("s_speedup", "mean")], yerr=agg[("s_speedup", "sem")],
                color=color, capsize=3, **version_style[v]
            )
            if v == "v4":
                last_x = threads.iloc[-1]
                last_y = agg[("s_speedup", "mean")].iloc[-1]
        if last_x is not None:
            label_points.append([org, color, last_x, last_y])
    if len(label_points) == 2:
        y0, y1 = label_points[0][3], label_points[1][3]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        min_gap = 0.06 * y_range
        gap = y1 - y0
        if abs(gap) < min_gap:
            shift = (min_gap - abs(gap)) / 2
            direction = 1 if gap >= 0 else -1
            label_points[0][3] -= shift * direction
            label_points[1][3] += shift * direction
    for org, color, x, y in label_points:
        ax.annotate(
            org, (x, y),
            textcoords="offset points", xytext=(8, 0),
            fontsize=9, color=color, va="center"
        )
    ax.axhline(1.0, color="grey", linewidth=0.8, linestyle=":")
    ax.set_title(mod.lower())
    ax.set_xlabel("threads")
    ax.set_ylabel("fold-speedup")
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin, xmax + 2)
    fig.canvas.draw()
    ticks = ax.get_xticks()
    labels = [t.get_text() for t in ax.get_xticklabels()]
    labels[-1] = ""
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
mem_ax = axes[len(modalities)]
sns.barplot(
    data=df, x="modality", y="max_rss", hue="version",
    ax=mem_ax, errorbar="se",
    palette={"v3": "tab:gray", "v4": "tab:red"}
)
mem_ax.set_title("memory usage")
mem_ax.set_xlabel("")
mem_ax.set_ylabel("max RSS (MB)")
mem_ax.set_xticklabels(mem_ax.get_xticklabels(), rotation=30, ha='right')
mem_ax.legend(fontsize=8, title=None, loc='upper left')
for j in range(len(modalities) + 1, len(axes)):
    axes[j].axis('off')

fig.savefig('benchmark.png', dpi=300)
fig.savefig('benchmark.pdf', dpi=300)
fig.savefig('benchmark.tiff', dpi=300)
