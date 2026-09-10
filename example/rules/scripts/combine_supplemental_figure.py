# Script to combine the 3 ATAC supplemental figures into a single figure

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Section 1: Read the images
fig1 = mpimg.imread('results/atac_pca.png')
fig2 = mpimg.imread('results/atac_correlation.png')
fig3 = mpimg.imread('results/atac_enrichment.png')


# Section 2: Create a new figure and add subplots
# Plot 1
ax1 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)
ax1.imshow(fig1)

# Plot 2
ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1)
ax2.imshow(fig2)

# Plot 3
ax3 = plt.subplot2grid((2, 2), (1, 1), colspan=1)
ax3.imshow(fig3)

# Add titles to each subplot
ax1.set_title('A', loc='left')
ax2.set_title('B', loc='left')
ax3.set_title('C', loc='left')

# remove axes for all subplots
for ax in [ax1, ax2, ax3]:
    ax.axis('off')

# save the combined figure
plt.savefig('results/supplemental_figure_combined.png', dpi=300, bbox_inches='tight')
plt.savefig('results/supplemental_figure_combined.pdf', dpi=300, bbox_inches='tight')
plt.savefig('results/supplemental_figure_combined.tiff', dpi=300, bbox_inches='tight')