rule bam_stats:
    input:
        bam = 'bamfiles/{cramfile}.bam'
    output:
        tsv = 'stats/{cramfile}.tsv'
    conda: 'env/sample_stats.yaml'
    threads: 2
    script:
        'scripts/sample_stats.py'

rule combine_stats:
    input:
        expand('stats/{cramfile}.tsv', cramfile=CRAMFILES)
    output:
        tsv = 'results/sample_stats.tsv'
    conda: 'env/sample_stats.yaml'
    localrule: True
    script:
        'scripts/combine_stats.py'
