rule multibamsummary_dt4:
  input:
    bam_files = lambda wildcards: [
        f"bamfiles/{i}.bam" for i in multibamSummary_samples[wildcards.run]
    ]
  output:
    npz = temp("output/mbs_{run}.dt4.npz")
  benchmark: repeat("benchmarks/multibamsummary/{run}_dt4.txt", config['repeats'])
  threads: 10
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
    npz = temp("output/mbs_{run}.dt3.npz")
  benchmark: repeat("benchmarks/multibamsummary/{run}_dt3.txt", config['repeats'])
  threads: 10
  resources:
    mem_mb = 20000,
    runtime = 1440
  shell:"""
  multiBamSummary_old bins -p {threads} \
    -o {output.npz} \
    -b {input.bam_files}
  """
