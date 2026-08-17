from pathlib import Path
import yaml
import os
import platform
import sys

if not config.get('os'):
    system_name = platform.system()
    if system_name == "Linux":
        release = platform.freedesktop_os_release()
        config['os'] = release.get("ID", "linux") or "linux"
    elif system_name == "Darwin":
        config['os'] = "macos"
    else:
        config['os'] = system_name.lower()

if not config.get('what'):
    config['what'] = 'all'

if not config.get('repeats'):
    config['repeats'] = 3
config['repeats'] = int(config['repeats'])

# Paths
repodir = Path(workflow.basedir)
with open(repodir / 'conf' / 'sources.yaml') as f:
    sampleconfig = yaml.safe_load(f)
BAMDIR = 'bamfiles'

CRAMFILES = [
    "human_chip_SRR28592124",
    "human_rna_SRR28012902",
    "human_wgs_SRR15494527",
    "triticum_chip_SRR1686799",
    "triticum_rna_SRR27822150",
    "triticum_wgs_SRR27887047",
]
CRAMFILESEXT = [f"{i}.bam" for i in CRAMFILES]
GENOMES = ["triticum", "human"]

wildcard_constraints:
    file = r".+\.(cram|fna\.gz|gtf\.gz|gtf|fna|bw)",
    n = r"\d+"

EXHAUSTIVE = bool(config.get('exhaustive', False))

def thread_range(max_threads):
    return list(range(1, max_threads + 1)) if EXHAUSTIVE else [max_threads]

ALLFILES = (
    expand("{f}.cram", f=CRAMFILES)
    + expand("{f}.bw", f=CRAMFILES)
    + expand("{g}.fna", g=GENOMES)
    + expand("{g}.gtf", g=GENOMES)
)

alignmentSieve_samples = {
    'alsieve_human_chip': 'human_chip_SRR28592124',
    'alsieve_human_rna': 'human_rna_SRR28012902',
    'alsieve_human_wgs': 'human_wgs_SRR15494527',
    'alsieve_triticum_chip': 'triticum_chip_SRR1686799',
    'alsieve_triticum_rna': 'triticum_rna_SRR27822150',
    'alsieve_triticum_wgs': 'triticum_wgs_SRR27887047'
}

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

ALSIEVE_THREADS = 4
BAMCOVERAGE_THREADS = 10
BAMCOMPARE_THREADS = 10
MULTIBAMSUMMARY_THREADS = 10
COMPUTEMATRIX_THREADS = 4

include: 'rules/alignmentsieve.smk'
include: 'rules/bamcoverage.smk'
include: 'rules/bamcompare.smk'
include: 'rules/multibamsummary.smk'
include: 'rules/computematrix.smk'
include: 'rules/plotter.smk'
include: 'rules/download_data.smk'
include: '../.shared_rules/sample_stats.smk'

def function_runners(conf_what):
    match conf_what:
        case 'alignmentSieve':
            return [
                expand("benchmarks/alignmentsieve/{alignmentsieve}_dt4_t{n}.txt", alignmentsieve=alignmentSieve_samples.keys(), n=thread_range(ALSIEVE_THREADS)),
                expand("benchmarks/alignmentsieve/{alignmentsieve}_dt3_t{n}.txt", alignmentsieve=alignmentSieve_samples.keys(), n=thread_range(ALSIEVE_THREADS)),
            ]
        case 'bamCoverage':
            return [
                expand("benchmarks/bamcoverage/{bamcoverage}_dt4_t{n}.txt", bamcoverage=bamCoverage_samples.keys(), n=thread_range(BAMCOVERAGE_THREADS)),
                expand("benchmarks/bamcoverage/{bamcoverage}_dt3_t{n}.txt", bamcoverage=bamCoverage_samples.keys(), n=thread_range(BAMCOVERAGE_THREADS)),
            ]
        case 'bamCompare':
            return [
                expand("benchmarks/bamcompare/{bamcompare}_dt4_t{n}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS)),
                expand("benchmarks/bamcompare/{bamcompare}_dt3_t{n}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS)),
            ]
        case 'multibamSummary':
            return [
                expand("benchmarks/multibamsummary/{multibamsummary}_dt4_t{n}.txt", multibamsummary=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS)),
                expand("benchmarks/multibamsummary/{multibamsummary}_dt3_t{n}.txt", multibamsummary=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS)),
            ]
        case 'computeMatrix':
            return [
                expand("benchmarks/computematrix/{computematrix}_dt4_t{n}.txt", computematrix=computeMatrix_samples.keys(), n=thread_range(COMPUTEMATRIX_THREADS)),
                expand("benchmarks/computematrix/{computematrix}_dt3_t{n}.txt", computematrix=computeMatrix_samples.keys(), n=thread_range(COMPUTEMATRIX_THREADS)),
            ]
        case _:
            return [
                expand("benchmarks/alignmentsieve/{alignmentsieve}_dt4_t{n}.txt", alignmentsieve=alignmentSieve_samples.keys(), n=thread_range(ALSIEVE_THREADS)),
                expand("benchmarks/alignmentsieve/{alignmentsieve}_dt3_t{n}.txt", alignmentsieve=alignmentSieve_samples.keys(), n=thread_range(ALSIEVE_THREADS)),
                expand("benchmarks/bamcoverage/{bamcoverage}_dt4_t{n}.txt", bamcoverage=bamCoverage_samples.keys(), n=thread_range(BAMCOVERAGE_THREADS)),
                expand("benchmarks/bamcoverage/{bamcoverage}_dt3_t{n}.txt", bamcoverage=bamCoverage_samples.keys(), n=thread_range(BAMCOVERAGE_THREADS)),
                expand("benchmarks/bamcompare/{bamcompare}_dt4_t{n}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS)),
                expand("benchmarks/bamcompare/{bamcompare}_dt3_t{n}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS)),
                expand("benchmarks/multibamsummary/{multibamsummary}_dt4_t{n}.txt", multibamsummary=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS)),
                expand("benchmarks/multibamsummary/{multibamsummary}_dt3_t{n}.txt", multibamsummary=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS)),
                'results/performance.csv',
                'results/performance.png'
            ]

rule all:
    input:
        expand("zenodo/{file}", file=ALLFILES),
        expand("bamfiles/{cramfile}.bam", cramfile=CRAMFILES),
        function_runners(config['what']),
        'results/sample_stats.tsv',
