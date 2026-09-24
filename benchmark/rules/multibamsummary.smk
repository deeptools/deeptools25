rule multibamsummary_dt4:
  input:
    bam_files = lambda wildcards: [
        f"bamfiles/{i}.bam" for i in multibamSummary_samples[wildcards.run]
    ]
  output:
    npz = temp("output/mbs_{run}.dt4_t{n}_rep{rep}.npz"),
    tsv = "output/mbs_{run}.dt4_t{n}_rep{rep}.tsv"
  benchmark: "benchmarks/multibamsummary/{run}_dt4_t{n}_rep{rep}.txt"
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 20000,
    runtime = 1440
  shell:"""
  multiBamSummary bins -p {threads} \
    -o {output.npz} \
    --outRawCounts {output.tsv} \
    -b {input.bam_files}
  """

rule multibamsummary_dt3:
  input:
    bam_files = lambda wildcards: [
        f"bamfiles/{i}.bam" for i in multibamSummary_samples[wildcards.run]
    ]
  output:
    npz = temp("output/mbs_{run}.dt3_t{n}_rep{rep}.npz"),
    tsv = "output/mbs_{run}.dt3_t{n}_rep{rep}.tsv"
  benchmark: "benchmarks/multibamsummary/{run}_dt3_t{n}_rep{rep}.txt"
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 20000,
    runtime = 1440
  shell:"""
  multiBamSummary_old bins -p {threads} \
    -o {output.npz} \
    --outRawCounts {output.tsv} \
    -b {input.bam_files} --genomeChunkSize 10000000
  """

rule compare_multibamsummary:
  input:
    dt3 = lambda wildcards: f"output/mbs_{wildcards.run}.dt3_t{MULTIBAMSUMMARY_THREADS}_rep1.tsv",
    dt4 = lambda wildcards: f"output/mbs_{wildcards.run}.dt4_t{MULTIBAMSUMMARY_THREADS}_rep1.tsv"
  output:
    tsv = "results/numerical_diff/multibamsummary_{run}.tsv"
  params:
    mode = "npz",
    tool = "multiBamSummary",
    sample = lambda wildcards: wildcards.run
  resources:
    mem_mb = 20000,
    runtime = 1440
  script:
    "scripts/compare_dt3_dt4.py"
