rule bamcoverage_dt4_human:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCoverage_samples[wildcards.bamcoverage]}.bam"
  output:
    bw = temp("output/bamcoverage_{bamcoverage}.dt4.bw")
  benchmark: repeat("benchmarks/bamcoverage/{bamcoverage}_dt4.txt", 3)
  threads: 10
  resources:
    mem_mb = 10000
  shell:'''
  bamCoverage -b {input.bam} -o {output.bw} \
    --binSize 10 --normalizeUsing RPKM -p {threads}
  '''

rule bamcoverage_dt3_human:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCoverage_samples[wildcards.bamcoverage]}.bam"
  output:
    bw = temp("output/bamcoverage_{bamcoverage}.dt3.bw")
  benchmark: repeat("benchmarks/bamcoverage/{bamcoverage}_dt3.txt", 3)
  threads: 10
  resources:
    mem_mb = 10000
  shell:'''
  bamCoverage_old -b {input.bam} -o {output.bw} \
    --binSize 10 --normalizeUsing RPKM -p {threads}
  '''
