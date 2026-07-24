from pathlib import Path
import yaml
import os
import platform
import sys

if not config.get('os'):
    if platform.system() == "Linux":
        config['os'] = platform.freedesktop_os_release().get("ID")
    else:
        config['os'] = platform.system().lower()

# Paths
repodir = Path(workflow.basedir)
with open(repodir / 'conf' / 'sources.yaml') as f:
    sampleconfig = yaml.safe_load(f)

CRAMFILES = [
    "human_chip_SRR28592124",
    "human_rna_SRR28012902",
    "human_wgs_SRR15494527",
    "triticum_chip_SRR1686799",
    "triticum_rna_SRR27822150",
    "triticum_wgs_SRR27887047",
]
GENOMES = ["triticum", "human"]

wildcard_constraints:
    file = r".+\.(cram|fna\.gz|gtf\.gz|gtf|fna|bw)"

ALLFILES = (
    expand("{f}.cram", f=CRAMFILES)
    + expand("{f}.bw", f=CRAMFILES)
    + expand("{g}.fna", g=GENOMES)
    + expand("{g}.gtf", g=GENOMES)
)

bamCoverage_samples = {
  'bcov_human_chip': 'human_chip_SRR28592124',
  'bcov_human_rna': 'human_rna_SRR28012902',
  'bcov_human_wgs': 'human_wgs_SRR15494527',
  'bcov_triticum_chip': 'triticum_chip_SRR1686799',
  'bcov_triticum_rna': 'triticum_rna_SRR27822150',
  'bcov_triticum_wgs': 'triticum_wgs_SRR27887047'
}

bamCompare_samples = {
  'bcom_human_chip': ['human_chip_SRR28592124', 'human_wgs_SRR15494527'],
  'bcom_triticum_chip': ['triticum_chip_SRR1686799', 'triticum_wgs_SRR27887047']
}

multibamSummary_samples = {
  'mbs_human_3': [
    'human_chip_SRR28592124',
    'human_rna_SRR28012902',
    'human_wgs_SRR15494527'
  ],
  'mbs_human_15': [
    'human_chip_SRR28592124',
    'human_rna_SRR28012902',
    'human_wgs_SRR15494527'
  ] * 5,
  'mbs_triticum_3': [
    'triticum_chip_SRR1686799',
    'triticum_rna_SRR27822150',
    'triticum_wgs_SRR27887047'
  ],
  'mbs_triticum_15': [
    'triticum_chip_SRR1686799',
    'triticum_rna_SRR27822150',
    'triticum_wgs_SRR27887047'
  ] * 5
}

computeMatrix_samples = {
  'cm_human_3': [
    'human_wgs_SRR15494527',
    'human_chip_SRR28592124',
    'human_rna_SRR28012902'
  ],
  'cm_human_15': [
    'human_wgs_SRR15494527',
    'human_chip_SRR28592124',
    'human_rna_SRR28012902'
  ] * 5,
  'cm_triticum_3': [
    'triticum_wgs_SRR27887047',
    'triticum_chip_SRR1686799',
    'triticum_rna_SRR27822150'
  ],
  'cm_triticum_15': [
    'triticum_wgs_SRR27887047',
    'triticum_chip_SRR1686799',
    'triticum_rna_SRR27822150'
  ] * 5
}

include: 'rules/bamcoverage.smk'
include: 'rules/bamcompare.smk'
include: 'rules/multibamsummary.smk'
include: 'rules/computematrix.smk'
include: 'rules/plotter.smk'
include: 'rules/download_data.smk'

rule all:
    input:
        expand("zenodo/{file}", file=ALLFILES),
        expand("bamfiles/{cramfile}.bam", cramfile=CRAMFILES),

        # bamCoverage
        expand("benchmarks/bamcoverage/{bamcoverage}_dt4.txt", bamcoverage=bamCoverage_samples.keys()),
        expand("benchmarks/bamcoverage/{bamcoverage}_dt3.txt", bamcoverage=bamCoverage_samples.keys()),
        # bamCompare
        expand("benchmarks/bamcompare/{bamcompare}_dt4.txt", bamcompare=bamCompare_samples.keys()),
        expand("benchmarks/bamcompare/{bamcompare}_dt3.txt", bamcompare=bamCompare_samples.keys()),
        # multiBamSummary
        expand("benchmarks/mbs/{run}_dt4.txt", run=multibamSummary_samples.keys()),
        expand("benchmarks/mbs/{run}_dt3.txt", run=multibamSummary_samples.keys()),
        # computeMatrix
        expand("benchmarks/computeMatrix/{run}_dt4.txt", run=computeMatrix_samples.keys()),
        expand("benchmarks/computeMatrix/{run}_dt3.txt", run=computeMatrix_samples.keys()),
        # performance
        'results/performance.csv',
        'results/performance.png'
