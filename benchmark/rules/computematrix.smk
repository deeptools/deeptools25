rule computeMatrix_dt4:
  input:
    bw_files = lambda wildcards: [
        f"zenodo/{i}.bw" for i in computeMatrix_samples[wildcards.run]
    ],
    gtf = lambda wildcards: f"zenodo/{wildcards.run.split('_')[1]}.gtf"
  output:
    npz = temp("output/computeMatrix_{run}.dt4_t{n}.npz")
  benchmark: repeat("benchmarks/computeMatrix/{run}_dt4_t{n}.txt", config['repeats'])
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 120000,
    runtime = 1440
  shell:'''
  computeMatrix reference-point -p {threads} \
    -o {output.npz} \
    -R {input.gtf} \
    -S {input.bw_files}
  '''

rule computeMatrix_dt3:
  input:
    bw_files = lambda wildcards: [
        f"zenodo/{i}.bw" for i in computeMatrix_samples[wildcards.run]
    ],
    gtf = lambda wildcards: f"zenodo/{wildcards.run.split('_')[1]}.gtf"
  output:
    npz = temp("output/computeMatrix_{run}.dt3_t{n}.npz")
  benchmark: repeat("benchmarks/computeMatrix/{run}_dt3_t{n}.txt", config['repeats'])
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 120000,
    runtime = 1440
  shell:'''
  computeMatrix_old reference-point -p {threads} \
    -o {output.npz} \
    -R {input.gtf} \
    -S {input.bw_files}
  '''
