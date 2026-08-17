rule multibamsummary_dt4:
  input:
    bam_files = lambda wildcards: [
        f"bamfiles/{i}.bam" for i in multibamSummary_samples[wildcards.run]
    ]
  output:
    npz = temp("output/mbs_{run}.dt4_t{n}.npz")
  benchmark: repeat("benchmarks/multibamsummary/{run}_dt4_t{n}.txt", config['repeats'])
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 20000,
    runtime = 1440
  shell:"""
  multiBamSummary bins -p {threads} \
    -o {output.npz} \
    -b {input.bam_files}
  """

rule multibamsummary_dt3:
  input:
    bam_files = lambda wildcards: [
        f"bamfiles/{i}.bam" for i in multibamSummary_samples[wildcards.run]
    ]
  output:
    npz = temp("output/mbs_{run}.dt3_t{n}.npz")
  benchmark: repeat("benchmarks/multibamsummary/{run}_dt3_t{n}.txt", config['repeats'])
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 20000,
    runtime = lambda wildcards: 2880 if int(wildcards.n) < 4 else 1440
  shell:"""
  multiBamSummary_old bins -p {threads} \
    -o {output.npz} \
    -b {input.bam_files}
  """
