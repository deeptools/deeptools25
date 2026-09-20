import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

df = pd.read_table(snakemake.input.csv)
df['organism'] = [i.split("_")[1] for i in df['sample']]
df['datatype'] = df['sample'].str.rsplit('_', n=1).str[1]
df['version'] = df['version'].map({'dt3': 'v3', 'dt4': 'v4'})
df['min'] = df['s'] / 60

modalities = df["modality"].unique()
organisms  = ["human", "triticum"]
versions   = ["v3", "v4"]

def sem(x):
    return np.std(x, ddof=1) / np.sqrt(len(x))

def summarize(sub, cols):
    return (
        sub.groupby("threads")[cols]
           .agg(['mean', sem])
           .reset_index()
    )

org_colors    = {"human": "#5B8FC7", "triticum": "#C76B95"}
version_style = {"v3": dict(linestyle="--", marker="s"), "v4": dict(linestyle="-",  marker="o")}

FIGSIZE = (12, 8)
RUNTIME_TICKS = [0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000]
MEMORY_TICKS = [50, 100, 1000, 5000, 10000, 50000]

def sparse_log_yaxis(ax, ticks):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())


def plot_org_lines(ax, sub_mod, ycol):
    label_points = []
    for org in organisms:
        sub_org = sub_mod[sub_mod["organism"] == org]
        color = org_colors[org]
        last = None
        for v in versions:
            sub = sub_org[sub_org["version"] == v]
            if sub.empty:
                continue
            agg = summarize(sub, [ycol])
            ax.errorbar(
                agg["threads"], agg[(ycol, "mean")], yerr=agg[(ycol, "sem")],
                color=color, capsize=3, **version_style[v]
            )
            if v == "v4":
                last = (agg["threads"].iloc[-1], agg[(ycol, "mean")].iloc[-1])
        if last is not None:
            label_points.append([org, color, *last, 0.0])
    annotate_labels(ax, label_points)


def annotate_labels(ax, label_points, min_gap_pt=11):
    if len(label_points) == 2:
        ax.autoscale_view()
        fig = ax.figure
        px = [ax.transData.transform((p[2], p[3]))[1] for p in label_points]
        gap_pt = (px[1] - px[0]) * 72 / fig.dpi
        if abs(gap_pt) < min_gap_pt:
            shift = (min_gap_pt - abs(gap_pt)) / 2
            direction = 1 if gap_pt >= 0 else -1
            label_points[0][4] = -shift * direction
            label_points[1][4] =  shift * direction
    for text, color, x, y, dy in label_points:
        ax.annotate(
            text, (x, y),
            textcoords="offset points", xytext=(8, dy),
            fontsize=9, color=color, va="center"
        )


def pad_xaxis_for_labels(ax, fig):
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin, xmax + 2)
    fig.canvas.draw()
    ticks = ax.get_xticks()
    labels = [t.get_text() for t in ax.get_xticklabels()]
    labels[-1] = ""
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)


def make_grid():
    n_panels = len(modalities) + 1
    ncols = -(-n_panels // 2)
    fig, axes = plt.subplots(
        nrows=2, ncols=ncols, figsize=FIGSIZE,
        squeeze=False, constrained_layout=True
    )
    return fig, axes.flatten()


def split_summary_panel(fig, slot, ycol, ylabel, title, ticks):
    """Replace `slot` by two stacked panels (human / triticum): open bars
    (organism colour, version line style) + individual samples, log y."""
    gs = slot.get_subplotspec().subgridspec(2, 1, hspace=0)
    slot.remove()
    ax_bot = fig.add_subplot(gs[1])
    ax_top = fig.add_subplot(gs[0], sharex=ax_bot)
    split_axes = dict(zip(organisms, [ax_top, ax_bot]))

    # one point per sample: mean over threads / replicates
    per_sample = (
        df.groupby(["modality", "version", "sample", "organism"], as_index=False)[ycol]
          .mean()
    )

    bar_width = 0.8
    hue_w = bar_width / len(versions)
    x_base = np.arange(len(modalities))
    rng = np.random.default_rng(0)

    for org, ax in split_axes.items():
        color = org_colors[org]
        sub_org = per_sample[per_sample.organism == org]
        # log axis: bars start from a common floor instead of 0
        floor = 10 ** np.floor(np.log10(sub_org[ycol].min() * 0.8))

        for v_i, v in enumerate(versions):
            xs = x_base - bar_width / 2 + hue_w * (v_i + 0.5)
            sub_v = sub_org[sub_org.version == v]
            stats = (sub_v.groupby("modality")[ycol]
                          .agg(["mean", sem]).reindex(modalities))
            ax.bar(xs, stats["mean"] - floor, bottom=floor, width=hue_w * 0.85,
                   fill=False, edgecolor=color, linewidth=1.2,
                   linestyle=version_style[v]["linestyle"], zorder=1)
            ax.errorbar(xs, stats["mean"], yerr=stats["sem"], fmt="none",
                        ecolor=color, capsize=3, linewidth=1, zorder=2)
            for m_i, mod in enumerate(modalities):
                pts = sub_v.loc[sub_v.modality == mod, ycol]
                jitter = rng.uniform(-hue_w * 0.22, hue_w * 0.22, len(pts))
                ax.scatter(xs[m_i] + jitter, pts, s=32, color=color,
                           marker=version_style[v]["marker"],
                           edgecolors="black", linewidths=0.6, alpha=0.9, zorder=3)

        sparse_log_yaxis(ax, ticks)
        ax.set_ylim(floor, sub_org[ycol].max() * 8)   # headroom for label / legend
        ax.set_ylabel(ylabel, fontsize=9)
        ax.text(0.02, 0.95, org, transform=ax.transAxes,
                color=color, fontsize=9, va="top")

    ax_top.set_title(title)
    ax_top.tick_params(labelbottom=False)
    ax_bot.set_xticks(x_base)
    ax_bot.set_xticklabels([m.lower() for m in modalities], rotation=30, ha="right")

    ver_handles = [
        Line2D([], [], color="black", label=v, linewidth=1.2, **version_style[v])
        for v in versions
    ]
    ax_top.legend(handles=ver_handles, fontsize=8, loc="upper right",
                  frameon=False, ncol=2, handlelength=2.5)


def line_figure(ycol, ylabel, ticks, summary_title, outname):
    fig, axes = make_grid()
    for i, mod in enumerate(modalities):
        ax = axes[i]
        sub_mod = df[df["modality"] == mod]
        ax.set_yscale("log")
        ax.set_ylim(sub_mod[ycol].min() / 1.5, sub_mod[ycol].max() * 1.5)
        plot_org_lines(ax, sub_mod, ycol)
        sparse_log_yaxis(ax, ticks)
        ax.set_title(mod.lower())
        ax.set_xlabel("threads")
        ax.set_ylabel(ylabel)
        pad_xaxis_for_labels(ax, fig)

    split_summary_panel(fig, axes[len(modalities)], ycol, ylabel, summary_title, ticks)

    for j in range(len(modalities) + 1, len(axes)):
        axes[j].axis('off')

    for ext in ("png", "pdf", "tiff"):
        fig.savefig(f'{outname}.{ext}', dpi=300)
    return fig

# Absolute runtimes per thread count
runfig = line_figure("min", "runtime (min)", RUNTIME_TICKS,"runtime (all threads)", "benchmark_runtime")
runfig.savefig(snakemake.output.runpng, dpi=300)
runfig.savefig(snakemake.output.runpdf, dpi=300)
runfig.savefig(snakemake.output.runtiff, dpi=300)

# Memory per thread count
memfig = line_figure("max_rss", "max RSS (MB)", MEMORY_TICKS,"memory usage", "benchmark_memory")
memfig.savefig(snakemake.output.mempng, dpi=300)
memfig.savefig(snakemake.output.mempdf, dpi=300)
memfig.savefig(snakemake.output.memtiff, dpi=300)
