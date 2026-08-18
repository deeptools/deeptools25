rule bam_stats:
    input:
        bam = BAMDIR + '/{cramfile}'
    output:
        tsv = 'stats/{cramfile}.tsv'
    conda: 'env/sample_stats.yml'
    threads: 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    script:
        'sample_stats.py'

rule combine_stats:
    input:
        expand('stats/{cramfile}.tsv', cramfile=CRAMFILESEXT)
    output:
        tsv = 'results/sample_stats.tsv'
    conda: 'env/sample_stats.yml'
    localrule: True
    script:
        'combine_stats.py'
