rule alsieve_dt4:
  input:
    bam = lambda wildcards: f"bamfiles/{alignmentSieve_samples[wildcards.alignmentsieve]}.bam"
  output:
    bam = temp("output/alignmentsieve_{alignmentsieve}.dt4.bam")
  benchmark: repeat("benchmarks/alignmentsieve/{alignmentsieve}_dt4.txt", 3)
  threads: 4
  resources:
    mem_mb = 10000
  shell:'''
  alignmentSieve -b {input.bam} -o {output.bam} \
    -p {threads} --ATACshift
  '''

rule alsieve_dt3:
  input:
    bam = lambda wildcards: f"bamfiles/{alignmentSieve_samples[wildcards.alignmentsieve]}.bam"
  output:
    bam = temp("output/alignmentsieve_{alignmentsieve}.dt3.bam")
  benchmark: repeat("benchmarks/alignmentsieve/{alignmentsieve}_dt3.txt", 3)
  threads: 4
  resources:
    mem_mb = 10000
  shell:'''
  alignmentSieve_old -b {input.bam} -o {output.bam} \
    -p {threads} --ATACshift
  '''
