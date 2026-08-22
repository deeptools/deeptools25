import string
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def add_panel(ax, path, label, title=None):
    img = mpimg.imread(path)
    ax.imshow(img)
    ax.axis('off')
    ax.text(0.0, 1.02, label, transform=ax.transAxes,
             fontsize=13, fontweight='bold', va='bottom', ha='left')
    if title:
        ax.set_title(title, fontsize=10, pad=12)


chips = snakemake.params.chips
n_chip = len(chips)
n_cols = max(n_chip, 3)

fig = plt.figure(figsize=(3 * n_cols, 9), constrained_layout=True)
gs = fig.add_gridspec(2, n_cols, height_ratios=[1, 1.3])
labels = iter(string.ascii_uppercase)

top_panels = [
    (snakemake.input.rna, 'RNA-seq DE'),
    (snakemake.input.atac, 'ATAC'),
    (snakemake.input.meth, 'CpG methylation'),
]
for col, (path, title) in enumerate(top_panels):
    ax = fig.add_subplot(gs[0, col])
    add_panel(ax, path, next(labels), title)
for col in range(len(top_panels), n_cols):
    fig.add_subplot(gs[0, col]).axis('off')

# bottom row: one panel per ChIP mark, in the order given by params.chips
for col, (chip, path) in enumerate(zip(chips, snakemake.input.chips)):
    ax = fig.add_subplot(gs[1, col])
    add_panel(ax, path, next(labels), chip)
for col in range(n_chip, n_cols):
    fig.add_subplot(gs[1, col]).axis('off')

fig.savefig(snakemake.output.pdf)
fig.savefig(snakemake.output.png, dpi=300)
