rule plot_benchmarks:
  localrule: True
  input:
    bcov4 = expand("benchmarks/bamcoverage/{bamcoverage}_dt4.txt", bamcoverage=bamCoverage_samples),
    bcov3 = expand("benchmarks/bamcoverage/{bamcoverage}_dt3.txt", bamcoverage=bamCoverage_samples),
    bcom4 = expand("benchmarks/bamcompare/{bamcompare}_dt4.txt", bamcompare=bamCompare_samples.keys()),
    bcom3 = expand("benchmarks/bamcompare/{bamcompare}_dt3.txt", bamcompare=bamCompare_samples.keys()),
    mbs4 = expand("benchmarks/mbs/{run}_dt4.txt", run=multibamSummary_samples.keys()),
    mbs3 = expand("benchmarks/mbs/{run}_dt3.txt", run=multibamSummary_samples.keys()),
    cm4 = expand("benchmarks/computeMatrix/{run}_dt4.txt", run=computeMatrix_samples.keys()),
    cm3 = expand("benchmarks/computeMatrix/{run}_dt3.txt", run=computeMatrix_samples.keys())
  output:
    csv = 'results/performance.csv',
    png = 'results/performance.png'
  params:
    os = config['os']
  conda: 'env/plotter.yaml'
  threads: 1
  script:
    'scripts/plotter.py'
