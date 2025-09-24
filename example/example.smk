from pathlib import Path
import yaml

# Paths
repodir = Path(workflow.basedir)
with open(repodir / 'conf' / 'example_sources.yaml') as f:
    sampleconfig = yaml.safe_load(f)
config['chromsizes'] = str(repodir / 'conf' / 'genome.chrom.sizes')
SAMPLES = sampleconfig['samples'].keys()
ATACSAMPLES = [sample for sample in SAMPLES if 'ATAC' in sample]
RNASAMPLES = [sample for sample in SAMPLES if 'RNA' in sample]
BSSAMPLES = [sample for sample in SAMPLES if 'BS' in sample]
CHIPS = list(set([sample.split('_')[3] for sample in sampleconfig['chipdict'].keys()]))
INH_CHIPS = {k: v for k, v in sampleconfig['chipdict'].items() if 'H3K27me3' in k or 'H3K9me3' in k}
INH_CHIP = ['H3K27me3', 'H3K9me3']

cmap = {
  'H3K4me3': 'Greens',
  'H3K27ac': 'Greens',
  'H3K4me1': 'Greens',
  'H3K9me3': 'Purples',
  'H3K27me3': 'Purples'
}

include: 'rules/get_data.smk'
include: 'rules/get_regions.smk'
include: 'rules/deeptools.smk'

rule all:
  input:
    # Get data in bam format. Road to this rule depends on zenodo or raw source.
    expand("deeptools_input/{sample}.bam", sample=SAMPLES),
    expand("deeptools_input/{sample}.bam.bai", sample=SAMPLES),
    expand("deeptools_input/{bssample}_CpG.bw", bssample=BSSAMPLES),
    'deeptools_input/mouse.fna',
    'deeptools_input/mouse.gtf',
    # Generate regions
    # 'deeptools_input/counts.txt',
    # 'deeptools_input/de_up.tsv',
    # 'deeptools_input/de_down.tsv',
    # expand('regions/{inh_chip}_uropa_finalhits.txt', inh_chip = INH_CHIP),
    'deeptools_input/upreg_tss.bed',
    'deeptools_input/downreg_tss.bed',
    'deeptools_input/upreg_genes.gtf',
    'deeptools_input/downreg_genes.gtf',
    'deeptools_input/upreg_H3K27me3.bed',
    'deeptools_input/downreg_H3K27me3.bed',
    'deeptools_input/upreg_H3K9me3.bed',
    'deeptools_input/downreg_H3K9me3.bed',
    # Deeptools rules
    # ChIP
    expand('deeptools_output/chip_{chip}.png', chip=CHIPS),
    'deeptools_output/atac.png',
    'deeptools_output/rna.png',
    'deeptools_output/meth.png',
    'deeptools_output/rna.png'
