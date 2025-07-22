from pathlib import Path
import yaml

# Paths
repodir = Path(workflow.basedir)
with open(repodir / 'conf' / 'example_sources.yaml') as f:
    sampleconfig = yaml.safe_load(f)
SAMPLES = sampleconfig['samples'].keys()

include: 'rules/get_data.smk'


rule all:
  input:
    # input data
     expand('fq/{sample}.valid', sample=SAMPLES),