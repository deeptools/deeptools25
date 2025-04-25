rule computeMatrix_dt4:
  params:
    bw_files = lambda wildcards: ' '.join(
        [f"zenodo/bigwigs/{i}.bw" for i in computeMatrix_samples[wildcards.run]]
    ),
    gtf = lambda wildcards: f"zenodo/gtf/{wildcards.run}.gtf"
  output:
    npz = temp("output/computeMatrix_{run}.dt4.npz")
  benchmark: repeat("benchmarks/computeMatrix/{run}_dt4.txt", 3)
  threads: 10
  shell:'''
  computeMatrix reference-point -p {threads} \
    -o {output.npz} \
    -R {params.gtf} \
    -S {params.bw_files}
  '''

rule computeMatrix_dt3:
  params:
    bw_files = lambda wildcards: ' '.join(
        [f"zenodo/bigwigs/{i}.bw" for i in computeMatrix_samples[wildcards.run]]
    ),
    gtf = lambda wildcards: f"zenodo/gtf/{wildcards.run}.gtf"
  output:
    npz = temp("output/computeMatrix_{run}.dt3.npz")
  benchmark: repeat("benchmarks/computeMatrix/{run}_dt3.txt", 3)
  threads: 10
  shell:'''
  computeMatrix_old reference-point -p {threads} \
    -o {output.npz} \
    -R {params.gtf} \
    -S {params.bw_files}
  '''