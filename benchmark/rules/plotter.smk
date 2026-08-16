rule plot_benchmarks:
  localrule: True
  input:
    alsieve4 = expand("benchmarks/alignmentsieve/{alignmentsieve}_dt4_t{n}.txt", alignmentsieve=alignmentSieve_samples, n=thread_range(ALSIEVE_THREADS)),
    alsieve3 = expand("benchmarks/alignmentsieve/{alignmentsieve}_dt3_t{n}.txt", alignmentsieve=alignmentSieve_samples, n=thread_range(ALSIEVE_THREADS)),
    bcov4 = expand("benchmarks/bamcoverage/{bamcoverage}_dt4_t{n}.txt", bamcoverage=bamCoverage_samples, n=thread_range(BAMCOVERAGE_THREADS)),
    bcov3 = expand("benchmarks/bamcoverage/{bamcoverage}_dt3_t{n}.txt", bamcoverage=bamCoverage_samples, n=thread_range(BAMCOVERAGE_THREADS)),
    bcom4 = expand("benchmarks/bamcompare/{bamcompare}_dt4_t{n}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS)),
    bcom3 = expand("benchmarks/bamcompare/{bamcompare}_dt3_t{n}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS)),
    mbs4 = expand("benchmarks/multibamsummary/{run}_dt4_t{n}.txt", run=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS)),
    mbs3 = expand("benchmarks/multibamsummary/{run}_dt3_t{n}.txt", run=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS)),
    cm4 = expand("benchmarks/computeMatrix/{run}_dt4_t{n}.txt", run=computeMatrix_samples.keys(), n=thread_range(COMPUTEMATRIX_THREADS)),
    cm3 = expand("benchmarks/computeMatrix/{run}_dt3_t{n}.txt", run=computeMatrix_samples.keys(), n=thread_range(COMPUTEMATRIX_THREADS))
  output:
    csv = 'results/performance.csv',
    png = 'results/performance.png'
  params:
    os = config['os']
  conda: 'env/plotter.yaml'
  threads: 1
  script:
    'scripts/plotter.py'
