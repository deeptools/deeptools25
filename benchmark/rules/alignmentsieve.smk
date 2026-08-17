rule alsieve_dt4:
  input:
    bam = lambda wildcards: f"bamfiles/{alignmentSieve_samples[wildcards.alignmentsieve]}.bam"
  output:
    bam = temp("output/alignmentsieve_{alignmentsieve}.dt4_t{n}_rep{rep}.bam")
  benchmark: "benchmarks/alignmentsieve/{alignmentsieve}_dt4_t{n}_rep{rep}.txt"
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 10000,
    runtime = 1440
  shell:'''
  alignmentSieve -b {input.bam} -o {output.bam} \
    -p {threads} --ATACshift
  '''

rule alsieve_dt3:
  input:
    bam = lambda wildcards: f"bamfiles/{alignmentSieve_samples[wildcards.alignmentsieve]}.bam"
  output:
    bam = temp("output/alignmentsieve_{alignmentsieve}.dt3_t{n}_rep{rep}.bam")
  benchmark: "benchmarks/alignmentsieve/{alignmentsieve}_dt3_t{n}_rep{rep}.txt"
  threads: lambda wildcards: int(wildcards.n)
  resources:
    mem_mb = 10000,
    runtime = 1440
  shell:'''
  alignmentSieve_old -b {input.bam} -o {output.bam} \
    -p {threads} --ATACshift
  '''
