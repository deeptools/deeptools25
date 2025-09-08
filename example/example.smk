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

include: 'rules/get_data.smk'
include: 'rules/deeptools.smk'

rule all:
  input:
    # Get data in bam format. Road to this rule depends on zenodo or raw source.
    expand("deeptools_input/{sample}.bam", sample=SAMPLES),
    expand("deeptools_input/{sample}.bam.bai", sample=SAMPLES),
    expand("deeptools_input/{bssample}_CpG.bedGraph", bssample=BSSAMPLES),
    # Deeptools rules
    # ChIP
    expand('deeptools_output/chip_{chip}.png', chip=CHIPS),
    'deeptools_output/atac.png',
    'deeptools_output/rna.png',
    #expand('deeptools_output/chip/{chipsample}.bw', chipsample=sampleconfig['chipdict'].keys()),
    #expand('deeptools_output/rna/{rnasample}.bw', rnasample=RNASAMPLES),
    #expand('deeptools_output/atac/{atacsample}.bw', atacsample=ATACSAMPLES),
