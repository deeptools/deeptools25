import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from deeptools import heatmapper
from matplotlib.colors import PowerNorm
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from matplotlib.path import Path
from matplotlib.patches import PathPatch

CHIPS = snakemake.params.chips
CMAP = snakemake.params.cmap
DE_COLORS = {'ns': '#e0e0e0', 'ns - viz': '#8c8c8c', 'down': '#2166ac', 'up': '#b2182b'}
GROUP_COLORS = {'up': DE_COLORS['up'], 'down': DE_COLORS['down'], 'non-de': DE_COLORS['ns - viz']}
GENE_COLORS = [DE_COLORS['up'], DE_COLORS['ns - viz'], DE_COLORS['down']]

######################### DE / HEATMAP defs #############################################
def parse_subreg(p):
    ixs = []
    with open(p) as f:
        for line in f:
            ix = line.strip().split('\t')[0]
            if ix != 'logFC':
                ixs.append(ix)
    return ixs


def plot_rna_de(ax, res_path, up_path, down_path, nonde_path):
    a = pd.read_table(res_path, sep='\t', index_col=0)

    up = parse_subreg(up_path)
    down = parse_subreg(down_path)
    nonde = parse_subreg(nonde_path)

    a['cat'] = 'ns'
    a.loc[up, 'cat'] = 'up'
    a.loc[down, 'cat'] = 'down'
    a.loc[nonde, 'cat'] = 'ns - viz'
    a['cat'] = pd.Categorical(a['cat'], categories=['ns', 'ns - viz', 'down', 'up'])

    zorders = {'ns': 1, 'ns - viz': 2, 'down': 3, 'up': 3}

    for cat in a['cat'].cat.categories:
        sub = a[a['cat'] == cat]
        ax.scatter(sub['logCPM'], sub['logFC'], c=DE_COLORS[cat], s=10, alpha=0.7,
                   zorder=zorders[cat], label=cat)

    ax.axhline(0, color='#e0e0e0', lw=0.8, linestyle='--', zorder=0)
    ax.set_xlabel('logCPM')
    ax.set_ylabel('logFC')

    n_up, n_down, n_nsviz = len(up), len(down), len(nonde)
    ax.text(0.98, 0.975, f'up: {n_up}', transform=ax.transAxes,
            ha='right', va='top', fontsize=10, color=DE_COLORS['up'])
    ax.text(0.98, 0.025, f'down: {n_down}', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=10, color=DE_COLORS['down'])
    ax.text(0.98, 0.55, f'ns - viz: {n_nsviz}', transform=ax.transAxes,
            ha='right', va='center', fontsize=10, color=DE_COLORS['ns - viz'])
    ax.set_xlim(-6, 20)

def load_matrix(path):
    hm = heatmapper.heatmapper()
    hm.read_matrix_file(path)
    return hm


def x_extent(hm):
    p = hm.parameters
    upstream = p['upstream'][0]
    downstream = p['downstream'][0]
    return -upstream, downstream

def nice_round(x, sig=2):
    if x == 0:
        return 0.0
    from math import log10, floor
    d = sig - int(floor(log10(abs(x)))) - 1
    return round(x, d)


def plot_mark_profiles_ridgeline(ax, marks, condition, height=0.85, show_ylabels=True):
    cond_key = 'ctrl' if condition == 'WT' else 'ko'

    for i, (mark, path) in enumerate(marks):
        hm = load_matrix(path)
        n_samples = hm.matrix.get_num_samples()
        n_groups = hm.matrix.get_num_groups()
        xmin, xmax = x_extent(hm)
        n_bins = hm.matrix.get_matrix(0, 0)['matrix'].shape[1]
        rel_pos = np.linspace(xmin, xmax, n_bins) / max(abs(xmin), abs(xmax))

        idxs = [s for s in range(n_samples) if cond_key in hm.matrix.sample_labels[s]]

        profiles = {}
        for g in range(n_groups):
            glabel = hm.matrix.group_labels[g]
            stack = np.stack([hm.matrix.get_matrix(g, s)['matrix'].filled(np.nan) for s in idxs], axis=0)
            avg_matrix = np.nanmean(stack, axis=0)
            profiles[glabel] = np.nanmean(avg_matrix, axis=0)

        all_vals = np.concatenate(list(profiles.values()))
        lo, hi = np.nanmin(all_vals), np.nanmax(all_vals)
        y0 = i

        for glabel, profile in profiles.items():
            if 'up' in glabel:
                color = GROUP_COLORS['up']
            elif 'down' in glabel:
                color = GROUP_COLORS['down']
            else:
                color = GROUP_COLORS['non-de']
            norm_profile = (profile - lo) / (hi - lo) if hi > lo else profile * 0
            ax.plot(rel_pos, y0 + norm_profile * height, color=color, lw=0.85)

    ax.axvline(0, color='#cccccc', lw=0.6, zorder=0)
    tick_positions = [i + height / 2 for i in range(len(marks))]
    ax.set_yticks(tick_positions, [m for m, _ in marks] if show_ylabels else [], fontsize=8)
    ax.tick_params(axis='y', left=True, labelleft=show_ylabels, length=3)
    ax.set_xlim(-1, 1)
    ax.set_xticks([])
    ax.text(0, -0.4, f'{condition}', ha='center', va='top')
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

def add_heatmap_panel(fig, outer_spec, matrix_path, mark_title, panel_label=None, cmap='YlOrRd',
                       zmin=None, zmax=None, pct=(5, 95), gamma=0.6,
                       missing_color='lightgrey', sample_titles=True,
                       sample_wspace=0, group_hspace=0.025, annot_width_ratio=0.12,
                       cbar_height_ratio=0.04, title_height_ratio=0.06):
    hm = load_matrix(matrix_path)
    n_samples = hm.matrix.get_num_samples()
    n_groups = hm.matrix.get_num_groups()
    xmin, xmax = x_extent(hm)

    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(missing_color)

    finite = np.ma.masked_invalid(hm.matrix.matrix).compressed()
    lo, hi = np.nanpercentile(finite, pct)
    zmin = lo if zmin is None else zmin
    zmax = hi if zmax is None else zmax
    # Give lower gamma to three narrow marks
    if mark_title in ['H3K27ac', 'H3K4me3', 'H3K4me1', 'H3K9me3', 'H3K27me3']:
        zmin = 0
    norm = PowerNorm(gamma=gamma, vmin=zmin, vmax=zmax)

    group_sizes = [hm.matrix.group_boundaries[g + 1] - hm.matrix.group_boundaries[g] for g in range(n_groups)]
    total = sum(group_sizes)
    title_size = total * title_height_ratio
    cbar_size = total * cbar_height_ratio
    width_ratios = [annot_width_ratio] + [1] * n_samples
    inner = outer_spec.subgridspec(n_groups + 2, n_samples + 1, wspace=sample_wspace, hspace=group_hspace,
                                    width_ratios=width_ratios,
                                    height_ratios=[title_size] + group_sizes + [cbar_size])

    ax_title = fig.add_subplot(inner[0, 1:])
    ax_title.axis('off')
    ax_title.set_title(mark_title)
    if panel_label:
        ax_title.text(-0.05, 1.06, panel_label, transform=ax_title.transAxes,
                       fontsize=13, fontweight='bold', va='bottom', ha='left')

    axes = []
    for g in range(n_groups):
        glabel = hm.matrix.group_labels[g]
        ax_annot = fig.add_subplot(inner[g + 1, 0])

        if 'up' in glabel:
            ax_annot.set_facecolor(GROUP_COLORS.get('up', '#cccccc'))
        elif 'down' in glabel:
            ax_annot.set_facecolor(GROUP_COLORS.get('down', '#cccccc'))
        else:
            ax_annot.set_facecolor(GROUP_COLORS.get('non-de', '#cccccc'))

        ax_annot.set_xticks([])
        ax_annot.set_yticks([])
        for spine in ax_annot.spines.values():
            spine.set_visible(False)

        for s in range(n_samples):
            ax = fig.add_subplot(inner[g + 1, s + 1])
            mat = hm.matrix.get_matrix(g, s)['matrix']
            ax.imshow(
                mat, aspect='auto', cmap=cm, norm=norm,
                interpolation='sinc', extent=(xmin, xmax, mat.shape[0], 0),
            )
            ax.set_yticks([])
            ax.set_xticks([])
            if sample_titles and g == 0:
                ax.set_title(hm.matrix.sample_labels[s].replace('-rep1.bw', '-1').replace('-rep2.bw', '-2').replace('rep', '').replace('ctrl', 'wt'), fontsize=6, ha='left', rotation=45, pad=1)
            axes.append(ax)

    sm = ScalarMappable(norm=norm, cmap=cm)
    ax_cbar = fig.add_subplot(inner[n_groups + 1, 1:])
    cbar = fig.colorbar(sm, cax=ax_cbar, orientation='horizontal')
    tick_pos = norm.inverse(np.linspace(0, 1, 5))
    cbar.set_ticks(tick_pos)
    cbar.set_ticklabels([f'{nice_round(v):g}' for v in tick_pos])
    cbar.ax.tick_params(labelsize=5, length=2, pad=1)

    return axes



################### Scheme defs ##################################################
def draw_track(ax, x0, x1, y_base, height, color, n=100, wiggle=0.3, seed=0, exon_regions=None, direction=1):
    x = np.linspace(x0, x1, n)
    rng = np.random.default_rng(seed)
    if exon_regions is not None:
        mask = np.zeros_like(x)
        for a, b in exon_regions:
            mask[(x >= a) & (x <= b)] = 1.0
        mask = np.convolve(mask, np.ones(7) / 7, mode='same')
    else:
        mask = np.sin(np.linspace(0, np.pi, n)) ** 0.5
    noise = wiggle * height * rng.standard_normal(n) * mask
    y = np.clip(height * mask + noise, 0.0, None)
    y = np.convolve(y, np.ones(5) / 5, mode='same')
    ax.fill_between(x, y_base, y_base + direction * y, color=color, alpha=0.7, lw=0)


def draw_wt_ko_tracks(ax, x0, x1, y, color, height=0.25, exon_regions=None, seed=0, labels=False):
    if labels:
        ax.text(x0 - 15, y + 0.16, 'WT', ha='right', va='bottom', fontsize=6, color='#555555')
        ax.text(x0 - 15, y - 0.16, 'KO', ha='right', va='top', fontsize=6, color='#555555')
    draw_track(ax, x0, x1, y + 0.15, height, color, exon_regions=exon_regions, seed=seed, direction=1)
    draw_track(ax, x0, x1, y - 0.15, height, color, exon_regions=exon_regions, seed=seed + 1, direction=-1)


def draw_gene(ax, y, color, name):
    exon1 = (-400, -100)
    intron = (-100, 100)
    exon2 = (100, 400)

    ax.annotate('', xy=(exon1[0], y + 0.3), xytext=(exon1[0], y - 0.1),
                arrowprops=dict(arrowstyle='-', color=color, lw=1.2, shrinkA=0, shrinkB=0))
    ax.annotate('', xy=(exon1[0] + 80, y + 0.3), xytext=(exon1[0], y + 0.3),
                arrowprops=dict(arrowstyle='-|>', color=color, lw=1.2, shrinkA=0, shrinkB=0))

    ax.add_patch(Rectangle((exon1[0], y - 0.1), exon1[1] - exon1[0], 0.2, color=color))
    ax.plot([intron[0], intron[1]], [y, y], color=color, lw=1.2)
    ax.add_patch(Rectangle((exon2[0], y - 0.1), exon2[1] - exon2[0], 0.2, color=color))

    ax.text(exon1[0], y - 0.15, name, ha='left', va='top', fontsize=8, color=color)


def make_regions(rng, x_start, x_end, n, min_gap=20, max_gap=150, min_w=30, max_w=160):
    regions = []
    x = x_start
    for _ in range(n):
        x += rng.integers(min_gap, max_gap)
        w = rng.integers(min_w, max_w)
        if x + w > x_end:
            break
        regions.append((x, w))
        x += w
    return regions


def draw_arc(ax, x0, x1, y0, bow, color, lw=1.5, alpha=0.8, label=None):
    xm = (x0 + x1) / 2
    ym = y0 + bow
    path = Path([(x0, y0), (xm, ym), (x1, y0)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    ax.add_patch(PathPatch(path, facecolor='none', edgecolor=color, lw=lw, alpha=alpha, zorder=1))
    if label:
        curve_ym = y0 + bow / 2
        va = 'bottom' if bow > 0 else 'top'
        ax.text(xm, curve_ym, label, ha='center', va=va, fontsize=6, color=color)


def draw_regions(ax, y, gene_color, seed=0, region_color='#bbbbbb', tss_x=-400, n_diff=1):
    rng = np.random.default_rng(seed)
    left_regions = make_regions(rng, -1050, -420, 3)
    right_regions = make_regions(rng, 420, 1050, 3)
    all_regions = left_regions + right_regions

    is_up = gene_color == DE_COLORS['up']
    is_down = gene_color == DE_COLORS['down']

    if is_up or is_down:
        candidates = list(range(len(all_regions)))
        n_sel = n_diff
    else:
        candidates = list(range(len(left_regions), len(all_regions)))  # right side only
        n_sel = 1

    sel_idxs = set(rng.choice(candidates, size=min(n_sel, len(candidates)), replace=False))

    for i, (x0, w) in enumerate(all_regions):
        is_sel = i in sel_idxs
        ax.add_patch(Rectangle((x0, y - 0.06), w, 0.12, color=region_color))

        if is_sel and is_up:
            wt_color, ko_color, wt_h, ko_h = '#bbbbbb', gene_color, 0.10, 0.30
        elif is_sel and is_down:
            wt_color, ko_color, wt_h, ko_h = gene_color, '#bbbbbb', 0.30, 0.10
        else:
            wt_color = ko_color = '#555555'
            wt_h = ko_h = 0.18

        draw_track(ax, x0, x0 + w, y + 0.15, wt_h, wt_color, seed=10 * (i + 1) + seed, direction=1)
        draw_track(ax, x0, x0 + w, y - 0.15, ko_h, ko_color, seed=10 * (i + 1) + seed + 1, direction=-1)

        if is_sel:
            if is_up:
                arc_y, bow, arc_color, lw, label = y + 0.35, 0.6, gene_color, 1.5, r'$max(log2FC_{ChIP})$'
            elif is_down:
                arc_y, bow, arc_color, lw, label = y - 0.3, -0.6, gene_color, 1.5, r'$max(log2FC_{ChIP})$'
            else:
                arc_y, bow, arc_color, lw, label = y - 0.3, -0.6, '#999999', 1.0, r'$min(log2FC_{ChIP})$'
            draw_arc(ax, tss_x, x0 + w / 2, arc_y, bow, arc_color, lw=lw, label=label)

    ax.text(-1000, y + 0.25, 'WT', ha='right', va='center', fontsize=6, color='#555555')
    ax.text(-1000, y - 0.25, 'KO', ha='right', va='center', fontsize=6, color='#555555')
    ax.text(-1000, y, 'regions', ha='right', va='center', fontsize=6, color='#555555')


################## FIGURE ##################################################

bottom_panels = ['atac', 'meth'] + CHIPS
n_cols = len(bottom_panels)

fig = plt.figure(figsize=(12, 8), tight_layout=True)
gs = fig.add_gridspec(2, 1, height_ratios=[0.8, 1])
gs_top = gs[0].subgridspec(1, 3, wspace=0.25, width_ratios=[1, 1.25, 1])
gs_bottom = gs[1].subgridspec(1, n_cols, wspace=0.15)

# MAPLOT ######################################################################################################
ax_rna = fig.add_subplot(gs_top[0, 0])
plot_rna_de(
    ax_rna,
    snakemake.input.rna_res,
    snakemake.input.rna_up,
    snakemake.input.rna_down,
    snakemake.input.rna_nonde,
)
ax_rna.set_title('RNA-seq DE', fontsize=10)

# SCHEME ######################################################################################################
ax_scheme = fig.add_subplot(gs_top[0, 1])
gene_names = ['gene A', 'gene B', 'gene C']
n_diffs = [1, 1, 1]
for y, color, name, nd in zip([2, 1, 0], GENE_COLORS, gene_names, n_diffs):
    draw_gene(ax_scheme, y, color, name)
    draw_regions(ax_scheme, y, color, seed=list(GENE_COLORS).index(color), n_diff=nd)

ax_scheme.set_xlim(-1150, 1100)
ax_scheme.set_ylim(-0.7, 2.7)
ax_scheme.set_xticks([])
ax_scheme.set_yticks([])
for spine in ax_scheme.spines.values():
    spine.set_visible(False)



# PROFILES ######################################################################################################
marks = [
    ('ATAC', snakemake.input.atac),
    ('Meth', snakemake.input.meth),
] + [(f.split('/')[-1].replace('.npz', '').replace('chip_', ''), f'{f}') for f in snakemake.input.chips]

gs_profiles_outer = gs_top[0, 2].subgridspec(2, 1, height_ratios=[0.08, 1], hspace=0.05)

ax_profiles_title = fig.add_subplot(gs_profiles_outer[0, 0])
ax_profiles_title.axis('off')
ax_profiles_title.set_title('Aggregate signal per group', fontsize=10)

gs_profiles = gs_profiles_outer[1, 0].subgridspec(1, 2, wspace=0.05)
ax_wt = fig.add_subplot(gs_profiles[0, 0])
ax_ko = fig.add_subplot(gs_profiles[0, 1])
plot_mark_profiles_ridgeline(ax_wt, marks, 'WT', show_ylabels=True)
plot_mark_profiles_ridgeline(ax_ko, marks, 'KO', show_ylabels=False)
ax_ko.set_ylim(ax_wt.get_ylim())

# HEATMAPS ######################################################################################################

add_heatmap_panel(fig, gs_bottom[0, 0], snakemake.input.atac, 'ATAC', cmap='Reds')
add_heatmap_panel(fig, gs_bottom[0, 1], snakemake.input.meth, 'CpG Meth.', cmap='Greys', missing_color='#ffffff')

for col, chip in enumerate(snakemake.input.chips, start=2):
    chipname = chip.split('/')[-1].replace('.npz', '').replace('chip_', '')
    add_heatmap_panel(fig, gs_bottom[0, col], chip, chipname, cmap=CMAP[chipname])

fig.canvas.draw()

y_top = max(ax_rna.get_position().y1, ax_scheme.get_position().y1, ax_profiles_title.get_position().y1) + 0.015
x_left = ax_rna.get_position().x0

fig.text(x_left - 0.02, y_top, 'A', fontsize=13, fontweight='bold', va='bottom', ha='left')
fig.text(ax_scheme.get_position().x0 - 0.02, y_top, 'B', fontsize=13, fontweight='bold', va='bottom', ha='left')
fig.text(ax_profiles_title.get_position().x0 - 0.02, y_top, 'C', fontsize=13, fontweight='bold', va='bottom', ha='left')
fig.text(x_left - 0.02, gs[1].get_position(fig).y1 + 0.01, 'D', fontsize=13, fontweight='bold', va='bottom', ha='left')

fig.savefig(snakemake.output.pdf)
fig.savefig(snakemake.output.png, dpi=300)
fig.savefig(snakemake.output.tiff, dpi=300)
