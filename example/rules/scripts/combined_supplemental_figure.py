# Script to combine the 3 ATAC supplemental figures into a single figure

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="compressed", figsize=(8, 8))
gs = GridSpec(2, 2, figure=fig, wspace=0, hspace=0)

# Section 1: Read the images
fig1 = mpimg.imread('results/chip_pca.png')
fig2 = mpimg.imread('results/chip_correlation.png')
fig3 = mpimg.imread('results/atac_enrichment.png')

# Section 2: Create a new figure and add subplots
# Plot 1
ax1 =  fig.add_subplot(gs[:, 0])
ax1.imshow(fig1)

# Plot 2
ax2 =  fig.add_subplot(gs[0, 1])
ax2.imshow(fig2)

# Plot 3
ax3 =  fig.add_subplot(gs[1, 1])
ax3.imshow(fig3)

# Add titles to each subplot
ax1.set_title('A', loc='left', size=8)
ax2.set_title('B', loc='left', size=8)
ax3.set_title('C', loc='left', size=8)

# remove axes for all subplots
for ax in [ax1, ax2, ax3]:
    ax.axis('off')

# save the combined figure
plt.savefig('results/supplemental_figure_combined.png', dpi=300, bbox_inches='tight')
plt.savefig('results/supplemental_figure_combined.pdf', dpi=300, bbox_inches='tight')
plt.savefig('results/supplemental_figure_combined.tiff', dpi=300, bbox_inches='tight')