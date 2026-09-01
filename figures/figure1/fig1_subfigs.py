from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT_DIR = Path(__file__).parent / "subfigs"

INK = "#111111"
RED = "#c0392b"
ORANGE = "#d98c3d"
BLUE = "#4a90d9"

GENOME_LEN = 100
BLACKLIST = (45, 62)
READ_H = 0.55
EXTRA_RED = 6
SAMPLE_GROUPS = [0, 0, 0, 1, 1, 1]
GROUP_COLORS = {0: RED, 1: BLUE}
ENRICHMENT_COLOR = "#9E4A06"
rng = np.random.default_rng(7)


def make_fragments(n=55, min_len=4, max_len=9, min_gap=4, max_gap=14, pair_prob=0.55):
    fragments = []
    for _ in range(n):
        paired = rng.random() < pair_prob
        r1_len = int(rng.integers(min_len, max_len))
        if paired:
            r2_len = int(rng.integers(min_len, max_len))
            gap = int(rng.integers(min_gap, max_gap))
            total = r1_len + gap + r2_len
            start = int(rng.integers(0, GENOME_LEN - total))
            reads = ((start, start + r1_len), (start + r1_len + gap, start + total))
            span = (start, start + total)
        else:
            start = int(rng.integers(0, GENOME_LEN - r1_len))
            reads = ((start, start + r1_len),)
            span = reads[0]
        fragments.append({"reads": reads, "span": span})

    fragments.sort(key=lambda f: f["span"][0])

    row_ends = []
    for frag in fragments:
        start, end = frag["span"]
        placed = False
        for i, row_end in enumerate(row_ends):
            if start > row_end + 2:
                frag["row"] = i
                row_ends[i] = end
                placed = True
                break
        if not placed:
            frag["row"] = len(row_ends)
            row_ends.append(end)

    return fragments

def classify(s, e):
    bl_start, bl_end = BLACKLIST
    if s >= bl_start and e <= bl_end:
        return "in"
    if s < bl_end and e > bl_start:
        return "partial"
    return "out"

def segments_for(r):
    s, e = r["s"], r["e"]
    if r["kind"] == "in" or r["extra_removed"]:
        return [(s, e, RED)]
    if r["kind"] == "partial":
        bl_start, bl_end = BLACKLIST
        segs = [(max(s, bl_start), min(e, bl_end), RED)]
        if s < bl_start:
            segs.append((s, min(e, bl_start), INK))
        if e > bl_end:
            segs.append((max(s, bl_end), e, INK))
        return segs
    return [(s, e, INK)]


def is_fully_removed(r):
    return r["kind"] == "in" or r["extra_removed"]

def build_reads(fragments):
    reads = []
    for frag_id, frag in enumerate(fragments):
        for s, e in frag["reads"]:
            reads.append({
                "frag": frag_id, "row": frag["row"], "s": s, "e": e,
                "kind": classify(s, e), "extra_removed": False,
            })

    candidates = [r for r in reads if r["kind"] == "out"]
    picks = rng.choice(len(candidates), size=min(EXTRA_RED, len(candidates)), replace=False)
    for i in picks:
        candidates[i]["extra_removed"] = True

    return reads


def draw_blacklist(ax, y0, y1):
    bl_start, bl_end = BLACKLIST
    ax.axvspan(bl_start, bl_end, color=RED, alpha=0.12, zorder=0, lw=0)
    for x in (bl_start, bl_end):
        ax.plot([x, x], [y0, y1], color=RED, lw=1.4, zorder=1)

def kept_segments(reads):
    segments = []
    for r in reads:
        if is_fully_removed(r):
            continue
        for s, e, color in segments_for(r):
            if color == INK:
                segments.append((s, e))
    return segments


def coverage_profile(segments):
    depth = np.zeros(GENOME_LEN, dtype=int)
    for s, e in segments:
        depth[s:e] += 1
    return depth


def smooth(depth, sigma=2.5):
    radius = int(sigma * 3)
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-x ** 2 / (2 * sigma ** 2))
    kernel /= kernel.sum()
    return np.convolve(depth, kernel, mode="same")


def plot(fragments, reads, path, keep_red):
    n_rows = max(f["row"] for f in fragments) + 1
    fig, ax = plt.subplots(figsize=(4.3, 2.9), dpi=300)
    draw_blacklist(ax, -0.7, n_rows - 0.3)

    by_frag = {}
    for r in reads:
        by_frag.setdefault(r["frag"], []).append(r)

    for frag_id, frag in enumerate(fragments):
        frag_reads = by_frag[frag_id]
        present = frag_reads if keep_red else [r for r in frag_reads if not is_fully_removed(r)]

        if len(frag_reads) == 2 and len(present) == 2:
            r1, r2 = sorted(present, key=lambda r: r["s"])
            ax.plot([r1["e"], r2["s"]], [frag["row"], frag["row"]],
                    color=INK, lw=0.9, ls="--", alpha=0.8, zorder=1)

        for r in present:
            for s, e, color in segments_for(r):
                if not keep_red and color == RED:
                    continue
                ax.add_patch(Rectangle(
                    (s, r["row"] - READ_H / 2), e - s, READ_H,
                    facecolor=color, edgecolor=color, zorder=2,
                ))

    ax.set_xlim(-2, GENOME_LEN + 2)
    ax.set_ylim(-0.9, n_rows - 0.1)
    ax.axis("off")
    fig.tight_layout(pad=0.15)
    fig.savefig(path, dpi=300, transparent=True)
    plt.close(fig)

def plot_bigwig(reads, path):
    depth = smooth(coverage_profile(kept_segments(reads)))
    x = np.arange(GENOME_LEN)
    y_top = depth.max() * 1.2

    fig, ax = plt.subplots(figsize=(4.3, 2.9), dpi=300)

    ax.fill_between(x, depth, color=ORANGE, alpha=0.35, zorder=1)
    ax.plot(x, depth, color=ORANGE, lw=1.3, zorder=2)

    ax.set_xlim(-2, GENOME_LEN + 2)
    ax.set_ylim(0, y_top)
    ax.axis("off")
    fig.tight_layout(pad=0.15)
    fig.savefig(path, dpi=300, transparent=True)
    plt.close(fig)

def plot_pca_scatter(path, seed=21):
    rng2 = np.random.default_rng(seed)
    n = len(SAMPLE_GROUPS)
    centers = [(-2.1, 0.6), (2.3, -0.5)]

    fig, ax = plt.subplots(figsize=(3.4, 3.5), dpi=300)

    for i in range(n):
        cx, cy = centers[SAMPLE_GROUPS[i]]
        x = cx + rng2.normal(0, 0.5)
        y = cy + rng2.normal(0, 0.5)
        ax.scatter(x, y, s=90, color=GROUP_COLORS[SAMPLE_GROUPS[i]], alpha=0.8, zorder=3)

    ax.set_xlabel("PC1 (46.2% of var. explained)", fontsize=7.5)
    ax.set_ylabel("PC2 (23.7% of var. explained)", fontsize=7.5)
    ax.set_title("PCA", fontsize=9)
    ax.tick_params(labelsize=7)
    fig.tight_layout(pad=0.4)
    fig.savefig(path, dpi=300, transparent=True)
    plt.close(fig)

def plot_scree(path):
    values = np.array([46.2, 23.7, 14.1, 8.4, 4.9, 2.7])
    ind = np.arange(1, len(values) + 1)

    fig, ax = plt.subplots(figsize=(3.4, 3.5), dpi=300)
    ax.plot(ind, values, "bo-", ms=4, lw=1, label="individual")
    ax.plot(ind, values.cumsum(), "ro-", ms=4, lw=1, label="accumulative")
    ax.axhline(1, color="#999999", lw=0.8, ls=":", zorder=1)
    ax.set_xticks(ind)
    ax.set_xlabel("Principal Component", fontsize=8)
    ax.set_ylabel("Variability", fontsize=8.5)
    ax.set_title("Scree plot", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6.5, labelcolor=["blue", "red"], frameon=False, loc="center right")
    fig.tight_layout(pad=0.4)
    fig.savefig(path, dpi=300, transparent=True)
    plt.close(fig)

def plot_corr_heatmap(path, seed=23):
    rng2 = np.random.default_rng(seed)
    n = len(SAMPLE_GROUPS)
    labels = [f"sample{i + 1}" for i in range(n)]

    corr = np.ones((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            same_group = SAMPLE_GROUPS[i] == SAMPLE_GROUPS[j]
            value = rng2.uniform(0.90, 0.98) if same_group else rng2.uniform(0.55, 0.72)
            corr[i, j] = value
            corr[j, i] = value

    fig = plt.figure(figsize=(4.0, 3.6), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[7, 1])
    ax_dendro = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_cbar = fig.add_subplot(gs[1, 1])
    dend_color = "darkred"

    def branch(x0, x1, y):
        ax_dendro.plot([x0, x1], [y, y], color=dend_color, lw=1)

    def joint(x, y0, y1):
        ax_dendro.plot([x, x], [y0, y1], color=dend_color, lw=1)

    leaf_y = [i + 0.5 for i in range(n)]
    a_y = (leaf_y[0] + leaf_y[1]) / 2
    b_y = (a_y + leaf_y[2]) / 2
    c_y = (leaf_y[3] + leaf_y[4]) / 2
    d_y = (c_y + leaf_y[5]) / 2

    branch(0, 0.3, leaf_y[0]); branch(0, 0.3, leaf_y[1]); joint(0.3, leaf_y[0], leaf_y[1])
    branch(0.3, 0.6, a_y); branch(0, 0.6, leaf_y[2]); joint(0.6, a_y, leaf_y[2])
    branch(0, 0.3, leaf_y[3]); branch(0, 0.3, leaf_y[4]); joint(0.3, leaf_y[3], leaf_y[4])
    branch(0.3, 0.6, c_y); branch(0, 0.6, leaf_y[5]); joint(0.6, c_y, leaf_y[5])
    branch(0.6, 1.0, b_y); branch(0.6, 1.0, d_y); joint(1.0, b_y, d_y)

    ax_dendro.set_xlim(1.05, -0.05)
    ax_dendro.set_ylim(n, 0)
    ax_dendro.axis("off")

    im = ax_heat.pcolormesh(corr, cmap="jet", vmin=0.5, vmax=1.0,
                              edgecolors="black", linewidth=0.4)
    ax_heat.invert_yaxis()
    ax_heat.set_xticks([i + 0.5 for i in range(n)])
    ax_heat.set_yticks([i + 0.5 for i in range(n)])
    ax_heat.xaxis.tick_top()
    ax_heat.set_xticklabels(labels, rotation=45, ha="left", fontsize=6.5)
    ax_heat.yaxis.tick_right()
    ax_heat.set_yticklabels(labels, fontsize=6.5)
    ax_heat.tick_params(length=0)

    fig.colorbar(im, cax=ax_cbar, orientation="horizontal")
    ax_cbar.tick_params(labelsize=6)

    fig.savefig(path, dpi=300, transparent=True)
    plt.close(fig)




def plot_enrichment_bar(path, seed=25):
    rng2 = np.random.default_rng(seed)
    features = ["peaks", "exons"]
    base_by_feature = {"peaks": 55, "exons": 30}
    n = len(SAMPLE_GROUPS)
    sample_labels = [f"sample{i + 1}" for i in range(n)]

    fig, axes = plt.subplots(1, len(features), figsize=(3.8, 3.3), dpi=300)
    for ax, feature in zip(axes, features):
        values = np.clip(base_by_feature[feature] + rng2.normal(0, 6, size=n), 2, 98)
        ax.bar(np.arange(n), values, width=0.5, color=ENRICHMENT_COLOR,
               edgecolor=ENRICHMENT_COLOR, alpha=0.9, zorder=3)
        ax.set_xticks(np.arange(n))
        ax.set_xticklabels(sample_labels, fontsize=6, rotation=90)
        ax.set_ylabel(f"% {feature}", fontsize=7.5)
        ax.set_ylim(0.0, 100.0)
        ax.tick_params(labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(pad=0.4, w_pad=1.4)
    fig.savefig(path, dpi=300, transparent=True)
    plt.close(fig)

def plot_heatmap_example(path, seed=27, n_regions=400, n_bins=80):
    rng2 = np.random.default_rng(seed)
    x = np.linspace(-1, 1, n_bins)
    base = np.exp(-(x ** 2) / (2 * 0.18 ** 2))
    scale = rng2.uniform(0.4, 1.3, size=(n_regions, 1))
    matrix = base[None, :] * scale + rng2.normal(0, 0.05, size=(n_regions, n_bins))
    matrix = np.clip(matrix, 0, None)
    matrix = matrix[matrix.sum(axis=1).argsort()[::-1]]
    profile = matrix.mean(axis=0)

    profile_color = plt.get_cmap("gnuplot")(0.0)

    fig, (ax_top, ax_main) = plt.subplots(
        2, 1, figsize=(3.7, 3.3), dpi=300,
        gridspec_kw={"height_ratios": [1, 3]},
    )

    ax_top.plot(x, profile, color=profile_color, lw=1.3, zorder=2)
    ax_top.fill_between(x, profile, color=profile_color, alpha=0.25, zorder=1)
    ax_top.set_xlim(-1, 1)
    ax_top.set_xticks([])
    ax_top.set_ylabel("signal", fontsize=6.5)
    ax_top.tick_params(labelsize=6)
    ax_top.spines[["top", "right"]].set_visible(False)

    im = ax_main.imshow(matrix, aspect="auto", cmap="RdYlBu",
                          extent=[-1, 1, n_regions, 0])
    ax_main.set_xticks([-1, 0, 1])
    ax_main.set_xticklabels(["-3.0", "TSS", "3.0Kb"], fontsize=6.5)
    ax_main.set_xlabel("gene distance (bp)", fontsize=6.5)
    ax_main.set_ylabel("regions", fontsize=7)
    ax_main.set_yticks([])
    cbar = fig.colorbar(im, ax=ax_main, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=6)

    fig.tight_layout(pad=0.4, h_pad=0.3)
    fig.savefig(path, dpi=300, transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    fragments = make_fragments()
    reads = build_reads(fragments)
    plot(fragments, reads, OUT_DIR / "bamfile.png", keep_red=True)
    plot(fragments, reads, OUT_DIR / "bamfile_sieve.png", keep_red=False)
    plot_bigwig(reads, OUT_DIR / "bigwig.png")
    plot_pca_scatter(OUT_DIR / "pca_scatter.png")
    plot_scree(OUT_DIR / "scree.png")
    plot_corr_heatmap(OUT_DIR / "corr_heatmap.png")
    plot_enrichment_bar(OUT_DIR / "enrichment_bar.png")
    plot_heatmap_example(OUT_DIR / "heatmap_example.png")
