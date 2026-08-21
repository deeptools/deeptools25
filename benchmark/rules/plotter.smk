rule plot_benchmarks:
  localrule: True
  input:
    alsieve4 = expand("benchmarks/alignmentsieve/{alignmentsieve}_dt4_t{n}_rep{rep}.txt", alignmentsieve=alignmentSieve_samples, n=thread_range(ALSIEVE_THREADS), rep=REPS),
    alsieve3 = expand("benchmarks/alignmentsieve/{alignmentsieve}_dt3_t{n}_rep{rep}.txt", alignmentsieve=alignmentSieve_samples, n=thread_range(ALSIEVE_THREADS), rep=REPS),
    bcov4 = expand("benchmarks/bamcoverage/{bamcoverage}_dt4_t{n}_rep{rep}.txt", bamcoverage=bamCoverage_samples, n=thread_range(BAMCOVERAGE_THREADS), rep=REPS),
    bcov3 = expand("benchmarks/bamcoverage/{bamcoverage}_dt3_t{n}_rep{rep}.txt", bamcoverage=bamCoverage_samples, n=thread_range(BAMCOVERAGE_THREADS), rep=REPS),
    bcom4 = expand("benchmarks/bamcompare/{bamcompare}_dt4_t{n}_rep{rep}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS), rep=REPS),
    bcom3 = expand("benchmarks/bamcompare/{bamcompare}_dt3_t{n}_rep{rep}.txt", bamcompare=bamCompare_samples.keys(), n=thread_range(BAMCOMPARE_THREADS), rep=REPS),
    mbs4 = expand("benchmarks/multibamsummary/{run}_dt4_t{n}_rep{rep}.txt", run=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS), rep=REPS),
    mbs3 = expand("benchmarks/multibamsummary/{run}_dt3_t{n}_rep{rep}.txt", run=multibamSummary_samples.keys(), n=thread_range(MULTIBAMSUMMARY_THREADS), rep=REPS),
    cm4 = expand("benchmarks/computeMatrix/{run}_dt4_t{n}_rep{rep}.txt", run=computeMatrix_samples.keys(), n=thread_range(COMPUTEMATRIX_THREADS), rep=REPS),
    cm3 = expand("benchmarks/computeMatrix/{run}_dt3_t{n}_rep{rep}.txt", run=computeMatrix_samples.keys(), n=thread_range(COMPUTEMATRIX_THREADS), rep=REPS)
  output:
    csv = 'results/performance.csv',
    png = 'results/performance.png',
    tiff = 'results/performance.tiff',
    pdf = 'results/performance.pdf'
  params:
    os = config['os']
  conda: 'env/plotter.yaml'
  threads: 2
  script:
    'scripts/plotter.py'
