rule bamcompare_dt4:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][0]}.bam",
    ctrlbam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][1]}.bam"
  output:
    bw = temp("output/bamcompare_{bamcompare}.dt4.bw")
  benchmark: repeat("benchmarks/bamcompare/{bamcompare}_dt4.txt", 3)
  resources:
    mem_mb = 20000
  threads: 10
  shell:'''
  bamCompare -b1 {input.bam} -b2 {input.ctrlbam} -o {output.bw} \
    --binSize 10 -p {threads}
  '''

rule bamcompare_dt3:
  input:
    bam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][0]}.bam",
    ctrlbam = lambda wildcards: f"bamfiles/{bamCompare_samples[wildcards.bamcompare][1]}.bam"
  output:
    bw = temp("output/bamcompare_{bamcompare}.dt3.bw")
  benchmark: repeat("benchmarks/bamcompare/{bamcompare}_dt3.txt", 3)
  threads: 10
  resources:
    mem_mb = 20000
  shell:'''
  bamCompare_old -b1 {input.bam} -b2 {input.ctrlbam} -o {output.bw} \
    --binSize 10 -p {threads}
  '''
