from pathlib import Path
import yaml

# Paths
repodir = Path(workflow.basedir)
with open(repodir / 'conf' / 'example_sources.yaml') as f:
    sampleconfig = yaml.safe_load(f)
SAMPLES = sampleconfig['samples'].keys()
ATACSAMPLES = [sample for sample in SAMPLES if 'ATAC' in sample]
RNASAMPLES = [sample for sample in SAMPLES if 'RNA' in sample]
BSSAMPLES = [sample for sample in SAMPLES if 'BS' in sample]

include: 'rules/get_data.smk'
include: 'rules/deeptools.smk'

rule all:
  input:
    # Get data in bam format. Road to this rule depends on zenodo or raw source.
    expand("deeptools_input/{sample}.bam", sample=SAMPLES),
    expand("deeptools_input/{sample}.bam.bai", sample=SAMPLES),
    expand("deeptools_input/{bssample}_CpG.bedGraph", bssample=BSSAMPLES),