import os
os.environ["TMPDIR"] = "/tmp"

bamCoverage_samples = {
  'bcov_benchmark_human_chip': 'benchmark_human_chip_SRR28592124',
  'bcov_benchmark_human_rna': 'benchmark_human_rna_SRR28012902',
  'bcov_benchmark_human_wgs': 'benchmark_human_wgs_SRR15494527',
  'bcov_benchmark_triticum_chip': 'benchmark_triticum_chip_SRR1686799',
  'bcov_benchmark_triticum_rna': 'benchmark_triticum_rna_SRR27822150',
  'bcov_benchmark_triticum_wgs': 'benchmark_triticum_wgs_SRR27887047'
}

bamCompare_samples = {
  'bcom_benchmark_human_chip': ['benchmark_human_chip_SRR28592124', 'benchmark_human_wgs_SRR15494527'],
  'bcom_benchmark_triticum_chip': ['benchmark_triticum_chip_SRR1686799', 'benchmark_triticum_wgs_SRR27887047']
}

multibamSummary_samples = {
  'mbs_benchmark_human_3': [
    'benchmark_human_chip_SRR28592124',
    'benchmark_human_rna_SRR28012902',
    'benchmark_human_wgs_SRR15494527'
  ],
  'mbs_benchmark_human_15': [
    'benchmark_human_chip_SRR28592124',
    'benchmark_human_rna_SRR28012902',
    'benchmark_human_wgs_SRR15494527'
  ] * 5,
  'mbs_benchmark_triticum_3': [
    'benchmark_triticum_chip_SRR1686799',
    'benchmark_triticum_rna_SRR27822150',
    'benchmark_triticum_wgs_SRR27887047'
  ],
  'mbs_benchmark_triticum_15': [
    'benchmark_triticum_chip_SRR1686799',
    'benchmark_triticum_rna_SRR27822150',
    'benchmark_triticum_wgs_SRR27887047'
  ] * 5
}

computeMatrix_samples = {
  'cm_benchmark_human_3': [
    'benchmark_human_wgs_SRR15494527',
    'benchmark_human_chip_SRR28592124',
    'benchmark_human_rna_SRR28012902'
  ],
  'cm_benchmark_human_15': [
    'benchmark_human_wgs_SRR15494527',
    'benchmark_human_chip_SRR28592124',
    'benchmark_human_rna_SRR28012902'
  ] * 5,
  'cm_benchmark_triticum_3': [
    'benchmark_triticum_wgs_SRR27887047',
    'benchmark_triticum_chip_SRR1686799',
    'benchmark_triticum_rna_SRR27822150'
  ],
  'cm_benchmark_triticum_15': [
    'benchmark_triticum_wgs_SRR27887047',
    'benchmark_triticum_chip_SRR1686799',
    'benchmark_triticum_rna_SRR27822150'
  ] * 5
}

include: 'rules/benchmark_bamcoverage.smk'
include: 'rules/benchmark_bamcompare.smk'
include: 'rules/benchmark_multibamsummary.smk'
include: 'rules/benchmark_computematrix.smk'
include: 'rules/plotter.smk'

rule all:
  input:
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
