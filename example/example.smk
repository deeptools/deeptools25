from pathlib import Path
import yaml

# Paths
repodir = Path(workflow.basedir)
with open(repodir / 'conf' / 'sources.yaml') as f:
    sampleconfig = yaml.safe_load(f)
config['chromsizes'] = str(repodir / 'conf' / 'genome.chrom.sizes')
config['rar'] = str(repodir / 'conf' / 'rar.bed')
with open(repodir / 'conf' / 'parameters.yaml') as f:
    config.update(yaml.safe_load(f))

# samples
ALLFILES = [
    'NPC_9sca_ctrl_H3K27ac_rep1.cram',
    'NPC_9sca_ctrl_H3K27ac_rep2.cram',
    'NPC_9sca_ctrl_H3K27me3_rep1.cram',
    'NPC_9sca_ctrl_H3K27me3_rep2.cram',
    'NPC_9sca_ctrl_H3K4me1_rep1.cram',
    'NPC_9sca_ctrl_H3K4me1_rep2.cram',
    'NPC_9sca_ctrl_H3K4me3_rep1.cram',
    'NPC_9sca_ctrl_H3K4me3_rep2.cram',
    'NPC_9sca_ctrl_H3K9me3_rep1.cram',
    'NPC_9sca_ctrl_H3K9me3_rep2.cram',
    'NPC_9sca_ctrl_input_rep1.cram',
    'NPC_9sca_ctrl_input_rep2.cram',
    'NPC_9sca_ctrl_RNA_rep1.cram',
    'NPC_9sca_ctrl_RNA_rep2.cram',
    'NPC_9sca_ctrl_ATAC_rep1.cram',
    'NPC_9sca_ctrl_ATAC_rep2.cram',
    'NPC_9sca_ctrl_BS_rep1.cram',
    'NPC_9sca_ctrl_BS_rep2.cram',
    'NPC_9sca_ko_ATAC_rep1.cram',
    'NPC_9sca_ko_ATAC_rep2.cram',
    'NPC_9sca_ko_BS_rep1.cram',
    'NPC_9sca_ko_BS_rep2.cram',
    'NPC_9sca_ko_H3K27ac_rep1.cram',
    'NPC_9sca_ko_H3K27ac_rep2.cram',
    'NPC_9sca_ko_H3K27me3_rep1.cram',
    'NPC_9sca_ko_H3K27me3_rep2.cram',
    'NPC_9sca_ko_H3K4me1_rep1.cram',
    'NPC_9sca_ko_H3K4me1_rep2.cram',
    'NPC_9sca_ko_H3K4me3_rep1.cram',
    'NPC_9sca_ko_H3K4me3_rep2.cram',
    'NPC_9sca_ko_H3K9me3_rep1.cram',
    'NPC_9sca_ko_H3K9me3_rep2.cram',
    'NPC_9sca_ko_input1_rep1.cram',
    'NPC_9sca_ko_input1_rep2.cram',
    'NPC_9sca_ko_RNA_rep1.cram',
    'NPC_9sca_ko_RNA_rep2.cram',
    'human.fna',
    'human.gtf',
    'mouse.fna',
    'mouse.gtf'
]
SAMPLES = [sample.replace('.cram', '') for sample in ALLFILES if 'fna' not in sample and 'gtf' not in sample]
ATACSAMPLES = [sample for sample in SAMPLES if 'ATAC' in sample]
RNASAMPLES = [sample for sample in SAMPLES if 'RNA' in sample]
BSSAMPLES = [sample for sample in SAMPLES if 'BS' in sample]
# ChIP samples, and types
CHIPS = list(set([sample.split('_')[3] for sample in sampleconfig['chipdict'].keys()]))
BROADMARKS = ['H3K27me3', 'H3K9me3']

cmap = {
  'H3K4me3': 'Greens',
  'H3K27ac': 'Greens',
  'H3K4me1': 'Greens',
  'H3K9me3': 'Purples',
  'H3K27me3': 'Purples'
}

wildcard_constraints:
    sample = "|".join(SAMPLES),

include: 'rules/get_data.smk'
include: 'rules/get_regions.smk'
include: 'rules/deeptools.smk'

rule all:
  input:
    # Download data
    expand("zenodo/{file}", file=ALLFILES),
    expand("deeptools_input/{sample}.bam", sample=SAMPLES),
    expand("deeptools_input/{sample}.bam.bai", sample=SAMPLES),
    expand("deeptools_input/{bssample}_CpG.bw", bssample=BSSAMPLES),
    'deeptools_input/mouse.fna',
    'deeptools_input/mouse.gtf',

    # Regions
    expand(
      'regions/{mergedpeak}_uropa_finalhits.txt',
      mergedpeak = ['ATAC'] + CHIPS
    ),

    # ChIP
    expand('deeptools_output/chip_{chip}.png', chip=CHIPS),
    'deeptools_output/atac.png',
    'deeptools_output/rna.png',
    'deeptools_output/meth.png',
    'deeptools_output/rna.png'
