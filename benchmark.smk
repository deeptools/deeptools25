
bamCoverage_samples = [
  'human_chip_SRR28592124',
  'human_rna_SRR28012902',
  'human_wgs_SRR15494527',
  'triticum_chip_SRR1686799',
  'triticum_rna_SRR27822150',
  'triticum_wgs_SRR27887047'
]

bamCompare_samples = {
  'human_chip_SRR28592124': 'human_wgs_SRR15494527',
  'triticum_chip_SRR1686799': 'triticum_wgs_SRR27887047'
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
  'human_3': [
    'human_wgs_SRR15494527',
    'human_chip_SRR28592124',
    'human_rna_SRR28012904'
  ],
  'human_15': [
    'human_wgs_SRR15494527',
    'human_chip_SRR28592124',
    'human_rna_SRR28012904'
  ] * 5,
  'triticum_3': [
    'triticum_wgs_mapq20',
    'triticum_chip_SRR1686799_nodup',
    'triticum_rna_SRR27822150_mapq10'
  ],
  'triticum_15': [
    'triticum_wgs_mapq20',
    'triticum_chip_SRR1686799_nodup',
    'triticum_rna_SRR27822150_mapq10'
  ] * 5
}

include: 'rules/benchmark_bamcoverage.smk'
include: 'rules/benchmark_bamcompare.smk'
include: 'rules/benchmark_multibamsummary.smk'
include: 'rules/benchmark_computematrix.smk'

rule all:
  input:
    # bamCoverage
    expand("benchmarks/bamcoverage/{bamcoverage}_dt4.txt", bamcoverage=bamCoverage_samples),
    expand("benchmarks/bamcoverage/{bamcoverage}_dt3.txt", bamcoverage=bamCoverage_samples),
    # bamCompare
    expand("benchmarks/bamcompare/{bamcompare}_dt4.txt", bamcompare=bamCompare_samples.keys()),
    expand("benchmarks/bamcompare/{bamcompare}_dt3.txt", bamcompare=bamCompare_samples.keys()),
    # multiBamSummary
    expand("benchmarks/mbs/{run}_dt4.txt", run=multibamSummary_samples.keys()),
    expand("benchmarks/mbs/{run}_dt3.txt", run=multibamSummary_samples.keys()),
    # computeMatrix
    expand("benchmarks/computeMatrix/{run}_dt4.txt", run=computeMatrix_samples.keys()),
    expand("benchmarks/computeMatrix/{run}_dt3.txt", run=computeMatrix_samples.keys())