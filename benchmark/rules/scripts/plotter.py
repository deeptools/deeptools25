import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, ScalarFormatter, NullFormatter
import re

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
versions   = ["v3", "v4"]

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

org_colors     = {"human": "#5B8FC7", "triticum": "#C76B95"}
version_style  = {"v3": dict(linestyle="--", marker="s"),
                  "v4": dict(linestyle="-",  marker="o")}
version_colors = {"v3": "tab:gray", "v4": "tab:red"}

MEM_LOG_Y = True

n_panels = len(modalities) + 1
ncols = -(-n_panels // 2)
fig, axes = plt.subplots(
    nrows=2, ncols=ncols, figsize=(12, 8),
    squeeze=False, constrained_layout=True
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
        for v in versions:
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

slot = axes[len(modalities)]
gs = slot.get_subplotspec().subgridspec(2, 1, hspace=0)
slot.remove()
ax_bot = fig.add_subplot(gs[1])
ax_top = fig.add_subplot(gs[0], sharex=ax_bot)
mem_axes = dict(zip(organisms, [ax_top, ax_bot]))

per_sample = (
    df.groupby(["modality", "version", "sample", "organism"], as_index=False)["max_rss"]
      .mean()
)

bar_width = 0.8
hue_w = bar_width / len(versions)
x_base = np.arange(len(modalities))
rng = np.random.default_rng(0)

for org, ax in mem_axes.items():
    color = org_colors[org]
    sub_org = per_sample[per_sample.organism == org]

    if MEM_LOG_Y:
        ymin_data = sub_org["max_rss"].min()
        floor = 10 ** np.floor(np.log10(ymin_data * 0.8))
    else:
        floor = 0

    for v_i, v in enumerate(versions):
        xs = x_base - bar_width / 2 + hue_w * (v_i + 0.5)
        sub_v = sub_org[sub_org.version == v]
        stats = (sub_v.groupby("modality")["max_rss"]
                      .agg(["mean", sem]).reindex(modalities))

        ax.bar(xs, stats["mean"] - floor, bottom=floor, width=hue_w * 0.85,
               fill=False, edgecolor=color, linewidth=1.2,
               linestyle=version_style[v]["linestyle"], zorder=1)
        ax.errorbar(xs, stats["mean"], yerr=stats["sem"], fmt="none",
                    ecolor=color, capsize=3, linewidth=1, zorder=2)
        for m_i, mod in enumerate(modalities):
            pts = sub_v.loc[sub_v.modality == mod, "max_rss"]
            jitter = rng.uniform(-hue_w * 0.22, hue_w * 0.22, len(pts))
            ax.scatter(xs[m_i] + jitter, pts, s=32, color=color,
                       marker=version_style[v]["marker"],
                       edgecolors="black", linewidths=0.6, alpha=0.9, zorder=3)

    if MEM_LOG_Y:
        ax.set_yscale("log")
        top = sub_org["max_rss"].max()
        ax.set_ylim(floor, top * 8)
        ax.yaxis.set_major_locator(FixedLocator([50, 100, 1000, 5000, 10000, 50000]))
        fmt = ScalarFormatter(); fmt.set_scientific(False)
        ax.yaxis.set_major_formatter(fmt)
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
    ax.set_ylabel("max RSS (MB)", fontsize=9)
    ax.text(0.02, 0.95, org, transform=ax.transAxes,
            color=color, fontsize=9, va="top")

ax_top.set_title("memory usage")
ax_top.tick_params(labelbottom=False)
ax_bot.set_xticks(x_base)
ax_bot.set_xticklabels([m.lower() for m in modalities], rotation=30, ha="right")

ver_handles = [
    Line2D([], [], color="black", label=v, linewidth=1.2, **version_style[v])
    for v in versions
]
ax_top.legend(handles=ver_handles, fontsize=8, loc="upper right",
              frameon=False, ncol=2, handlelength=2.5)

for j in range(len(modalities) + 1, len(axes)):
    axes[j].axis('off')

fig.savefig(snakemake.output.png, dpi=300)
fig.savefig(snakemake.output.pdf, dpi=300)
fig.savefig(snakemake.output.tiff, dpi=300)
