rule bamcompare_dt4:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][0]}.bam",
    ctrlbam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][1]}.bam"
  output:
    bw = temp("output/bamcompare_{bamcompare}.dt4_t{n}.bw")
  benchmark: repeat("benchmarks/bamcompare/{bamcompare}_dt4_t{n}.txt", config['repeats'])
  resources:
    mem_mb = 20000
  threads: lambda wildcards: int(wildcards.n)
  shell:'''
  bamCompare -b1 {input.bam} -b2 {input.ctrlbam} -o {output.bw} \
    --binSize 10 -p {threads}
  '''

rule bamcompare_dt3:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][0]}.bam",
    ctrlbam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][1]}.bam"
  output:
    bw = temp("output/bamcompare_{bamcompare}.dt3_t{n}.bw")
  benchmark: repeat("benchmarks/bamcompare/{bamcompare}_dt3_t{n}.txt", config['repeats'])
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 20000
  shell:'''
  bamCompare_old -b1 {input.bam} -b2 {input.ctrlbam} -o {output.bw} \
    --binSize 10 -p {threads}
  '''
