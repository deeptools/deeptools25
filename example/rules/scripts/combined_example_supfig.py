import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from matplotlib.lines import Line2D
from deeptools.correlation import Correlation
from deeptools.sumCoveragePerBin import SumCoveragePerBin
from deeptools.utilities import smartLabels

NPZ = snakemake.input.npz
BAMS = snakemake.input.bams

MARK_COLORS = {"H3K27ac": "C1", "H3K27me3": "C2", "H3K4me1": "C3", "H3K4me3": "C4", "H3K9me3": "C5"}


def load_labels(npz_path):
    raw = np.load(npz_path)["labels"]
    if raw.dtype == np.uint8:
        return [bytes(row).decode("utf-8").rstrip("\x00") for row in raw]
    return [str(x) for x in raw]


def color_for(label):
    for mark, color in MARK_COLORS.items():
        if mark in label:
            return color
    return "C0"


def marker_for(label):
    return "o" if "ctrl" in label else "s"


def linestyle_for(label):
    return "-" if "ctrl" in label else "--"


def heatmap_label(label):
    """'NPC_9sca_ctrl_H3K27ac_rep1.bam' -> 'H3K27ac_wt-1'."""
    stem = label.replace(".bam", "")
    parts = stem.split("_")
    mark, rep = parts[-2], parts[-1]
    rep_num = rep.replace("rep", "")
    condition = "wt" if "ctrl" in label else "ko"
    return f"{mark}_{condition}-{rep_num}"


fig = plt.figure(figsize=(12, 10), tight_layout=True)
gs = fig.add_gridspec(2, 2, wspace=0.3, hspace=0.3)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axD = fig.add_subplot(gs[1, 1])

npz_labels = load_labels(NPZ)
colors = [color_for(l) for l in npz_labels]
markers = [marker_for(l) for l in npz_labels]

corr = Correlation(NPZ, labels=npz_labels)
corr.transpose = True
corr.rowCenter = False
corr.log2 = True
corr.ntop = 20000

Wt, eigenvalues = corr.plot_pca(plot_filename=None, PCs=[1, 2], cols=colors, marks=markers)
pvar = eigenvalues / eigenvalues.sum()

for i, label in enumerate(npz_labels):
    axA.scatter(Wt[0, i], Wt[1, i], color=colors[i], marker=markers[i], alpha=0.8, s=40, zorder=i)
axA.set_xlabel(f"PC1 ({100 * pvar[0]:.1f}% of var. explained)")
axA.set_ylabel(f"PC2 ({100 * pvar[1]:.1f}% of var. explained)")
axA.set_title("PCA", fontsize=10)
axA.spines[["top", "right"]].set_visible(False)

mark_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=color,
                      markersize=7, label=mark) for mark, color in MARK_COLORS.items()]
mark_legend = axA.legend(handles=mark_handles, fontsize=7, frameon=False, loc="upper left")
axA.add_artist(mark_legend)

geno_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="black", markersize=7, label="WT"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="black", markersize=7, label="KO"),
]
axA.legend(handles=geno_handles, fontsize=7, frameon=False, loc="lower right")
axA.margins(0.15)

n_scree = min(10, len(pvar))
ind = np.arange(1, n_scree + 1)
axB.plot(ind, pvar[:n_scree], "o-", color="blue", label="individual")
axB.plot(ind, pvar.cumsum()[:n_scree], "o-", color="red", label="accumulative")
axB.set_xticks(ind)
axB.set_xlabel("Principal component")
axB.set_ylabel("Variability")
axB.set_title("Scree plot", fontsize=10)
axB.legend(fontsize=8, frameon=False)
axB.spines[["top", "right"]].set_visible(False)

corr2 = Correlation(NPZ, corr_method="pearson", labels=npz_labels)
corr_matrix = corr2.compute_correlation()
labels2 = [heatmap_label(l) for l in corr2.labels]

gs_c = gs[1, 0].subgridspec(1, 2, width_ratios=[0.08, 1], wspace=0.02)
ax_dendro = fig.add_subplot(gs_c[0, 0])
ax_heat = fig.add_subplot(gs_c[0, 1])

linkage = sch.linkage(corr_matrix, method="centroid")
order = sch.dendrogram(linkage, orientation="left", no_plot=True)["leaves"]
sch.dendrogram(linkage, orientation="left", ax=ax_dendro, link_color_func=lambda k: "#888888")
ax_dendro.axis("off")

ordered = corr_matrix[np.ix_(order, order)]
im = ax_heat.pcolormesh(ordered, cmap="RdBu_r", vmin=0, vmax=1)
ax_heat.set_xlim(0, len(labels2))
ax_heat.set_ylim(0, len(labels2))
ax_heat.set_xticks(np.arange(len(order)) + 0.5)
ax_heat.set_xticklabels(np.array(labels2)[order], rotation=90, fontsize=6)
ax_heat.set_yticks([])
fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.02).ax.tick_params(labelsize=6)
ax_heat.set_title("Sample correlation (Pearson)", fontsize=10)

bam_labels = smartLabels(BAMS)

cr = SumCoveragePerBin(
    BAMS,
    binLength=500,
    numberOfSamples=500_000,
    numberOfProcessors=20,
)
num_reads_per_bin = cr.run()

total = num_reads_per_bin.shape[0]
x = np.arange(total, dtype=float) / total

for i, reads in enumerate(num_reads_per_bin.T):
    count = np.cumsum(np.sort(reads))
    count = count / count[-1]
    axD.plot(x, count, color=color_for(bam_labels[i]), linestyle=linestyle_for(bam_labels[i]))

axD.set_xlabel("rank")
axD.set_ylabel("fraction w.r.t. bin with highest coverage")
axD.set_title("Fingerprint", fontsize=10)
legend_handles = [
    Line2D([0], [0], color="black", linestyle="-", label="WT"),
    Line2D([0], [0], color="black", linestyle="--", label="KO"),
]
axD.legend(handles=legend_handles, fontsize=8, frameon=False)
axD.spines[["top", "right"]].set_visible(False)

fig.canvas.draw()
y_top = max(axA.get_position().y1, axB.get_position().y1) + 0.015
fig.text(axA.get_position().x0 - 0.02, y_top, "A", fontsize=13, fontweight="bold", va="bottom", ha="left")
fig.text(axB.get_position().x0 - 0.02, y_top, "B", fontsize=13, fontweight="bold", va="bottom", ha="left")
fig.text(ax_dendro.get_position().x0 - 0.02, ax_dendro.get_position().y1 + 0.015, "C",
        fontsize=13, fontweight="bold", va="bottom", ha="left")
fig.text(axD.get_position().x0 - 0.02, axD.get_position().y1 + 0.015, "D",
        fontsize=13, fontweight="bold", va="bottom", ha="left")

fig.savefig(snakemake.output.pdf)
fig.savefig(snakemake.output.png, dpi=300)
fig.savefig(snakemake.output.tiff, dpi=300)
