rule bamcompare_chip:
    input:
        bam = 'deeptools_input/{chipsample}.bam'
    output:
        bw = 'deeptools_output/{chipsample}.log2.bw'
    threads: 10
    params:
        input = lambda wildcards: f'deeptools_input/{sampleconfig['chipdict'][wildcards.chipsample]}.bam'
    shell:'''
    bamCompare -b1 {input.bam} -b2 {params.input} \
      -bs 10 -p {threads} -o {output.bw}
    '''

rule bamCoverage_atac:
    input:
        bam = 'deeptools_input/{atacsample}.bam'
    output:
        bw = 'deeptools_output/{atacsample}.bw'
    threads: 10
    shell:'''
    bamCoverage -b {input.bam} -o {output.bw} \
      --normalizeUsing RPKM -bs 10 -p {threads}
    '''

rule bamCoverage_RNA:
    input:
        bam = 'deeptools_input/{rnasample}.bam'
    output:
        bw = 'deeptools_output/{rnasample}.bw'
    threads: 10
    shell:'''
    bamCoverage -b {input.bam} -o {output.bw} \
      --normalizeUsing RPKM -bs 10 -p {threads}
    '''

rule bigwig_BS:
    input:
        bam = 'deeptools_input/{bssample}.bam'
    output:
        bw = 'deeptools_output/{bssample}.bw'
