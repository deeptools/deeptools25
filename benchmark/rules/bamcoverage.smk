rule bamcoverage_dt4_human:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCoverage_samples[wildcards.bamcoverage]}.bam"
  output:
    bw = temp("output/bamcoverage_{bamcoverage}.dt4_t{n}.bw")
  benchmark: repeat("benchmarks/bamcoverage/{bamcoverage}_dt4_t{n}.txt", config['repeats'])
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 16000,
    runtime = 1440
  shell:'''
  bamCoverage -b {input.bam} -o {output.bw} \
    --binSize 10 --normalizeUsing RPKM -p {threads}
  '''

rule bamcoverage_dt3_human:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCoverage_samples[wildcards.bamcoverage]}.bam"
  output:
    bw = temp("output/bamcoverage_{bamcoverage}.dt3_t{n}.bw")
  benchmark: repeat("benchmarks/bamcoverage/{bamcoverage}_dt3_t{n}.txt", config['repeats'])
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 10000,
    runtime = 1440
  shell:'''
  bamCoverage_old -b {input.bam} -o {output.bw} \
    --binSize 10 --normalizeUsing RPKM -p {threads}
  '''
