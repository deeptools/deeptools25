rule bam_stats:
    input:
        bam = 'deeptools_input/{sample}.bam',
        bai = 'deeptools_input/{sample}.bam.bai'
    output:
        tsv = 'stats/{sample}.tsv'
    conda: 'env/sample_stats.yaml'
    threads: 2
    script:
        'scripts/sample_stats.py'

rule combine_stats:
    input:
        expand('stats/{sample}.tsv', sample=SAMPLES)
    output:
        tsv = 'results/sample_stats.tsv'
    conda: 'env/sample_stats.yaml'
    localrule: True
    script:
        'scripts/combine_stats.py'
