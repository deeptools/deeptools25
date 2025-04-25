rule multibamsummary_dt4:
  params:
    bam_files = lambda wildcards: ' '.join(
        [f"zenodo/{i}.bam" for i in multibamSummary_samples[wildcards.run]]
    )
  output:
    npz = temp("output/mbs_{run}.dt4.npz")
  benchmark: repeat("benchmarks/mbs/{run}_dt4.txt", 3)
  threads: 10
  shell:'''
  multiBamSummary bins -p {threads} \
    -o {output.npz} \
    -b {params.bam_files}
  '''

rule multibamsummary_dt3:
  params:
    bam_files = lambda wildcards: ' '.join(
        [f"zenodo/{i}.bam" for i in multibamSummary_samples[wildcards.run]]
    )
  output:
    npz = temp("output/mbs_{run}.dt3.npz")
  benchmark: repeat("benchmarks/mbs/{run}_dt3.txt", 3)
  threads: 10
  shell:'''
  multiBamSummary_old bins -p {threads} \
    -o {output.npz} \
    -b {params.bam_files}
  '''
