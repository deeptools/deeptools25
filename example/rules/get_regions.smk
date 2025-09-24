rule get_counts:
    input:
        bams = expand('deeptools_input/{sample}.bam', sample=RNASAMPLES),
        bais = expand('deeptools_input/{sample}.bam.bai', sample=RNASAMPLES),
        gtf = 'deeptools_input/mouse.gtf'
    output:
        counts = 'regions/counts.txt'
    threads: 10
    shell:'''
    featureCounts -T {threads} -p -o {output.counts} \
      -a {input.gtf} {input.bams}
    '''

rule DE:
    input:
        counts = 'regions/counts.txt'
    output:
        up = 'regions/de_up.tsv',
        down = 'regions/de_down.tsv',
    params:
        padj = 0.1,
        l2fc = 2
    script:
        'scripts/DESeq2.R'

rule call_peaks:
    input:
        bam = 'deeptools_input/{sample}.bam',
        ctrl = lambda wildcards: f'deeptools_input/{INH_CHIPS[wildcards.sample]}.bam'
    output:
        peak = temp('regions/{sample}_peaks.broadPeak'),
        gappedpeak = temp('regions/{sample}_peaks.gappedPeak'),
        xls = temp('regions/{sample}_peaks.xls'),
    shell:'''
    macs3 callpeak --broad -q 1e-2 -t {input.bam} -c {input.ctrl} \
      --keep-dup all \
      --outdir regions \
      -n {wildcards.sample} -f BAMPE -g mm
    '''

rule merge_peaks:
    input:
        expand('regions/{sample}_peaks.broadPeak', sample=INH_CHIPS.keys())
    output:
        bed = 'regions/{inh_chip}.bed',
    params:
        peaks = lambda wildcards: expand('regions/{sample}_peaks.broadPeak', sample=[s for s in INH_CHIPS.keys() if wildcards.inh_chip in s]),
    shell:'''
    cat {params.peaks} | sort -k1,1 -k2,2n | bedtools merge > {output.bed}
    '''

rule annotate_peaks:
    input:
        gtf = 'deeptools_input/mouse.gtf',
        bed = 'regions/{inh_chip}.bed',
    output:
        beda = temp('regions/{inh_chip}_uropa_allhits.bed'),
        txta = temp('regions/{inh_chip}_uropa_allhits.txt'),
        bedf = temp('regions/{inh_chip}_uropa_finalhits.bed'),
        txtf = 'regions/{inh_chip}_uropa_finalhits.txt',
        json = temp('regions/{inh_chip}_uropa.json'),
        pdf = temp('regions/{inh_chip}_uropa_summary.pdf')
    threads: 10
    shell:'''
    uropa -b {input.bed} -g {input.gtf} --summary \
      --feature gene --distance 100000 100000 \
      --internals 1 -p {wildcards.inh_chip}_uropa -o regions \
      --show-attributes gene_id gene_name
    '''

rule parse_regions:
    input:
        up = 'regions/de_up.tsv',
        down = 'regions/de_down.tsv',
        gtf = 'deeptools_input/mouse.gtf',
        uro_k27 = 'regions/H3K27me3_uropa_finalhits.txt',
        uro_k9 = 'regions/H3K9me3_uropa_finalhits.txt'
    output:
        upbed = 'deeptools_input/upreg_tss.bed',
        downbed = 'deeptools_input/downreg_tss.bed',
        upgtf = 'deeptools_input/upreg_genes.gtf',
        downgtf = 'deeptools_input/downreg_genes.gtf',
        k27_up = 'deeptools_input/upreg_H3K27me3.bed',
        k27_down = 'deeptools_input/downreg_H3K27me3.bed',
        k9_up = 'deeptools_input/upreg_H3K9me3.bed',
        k9_down = 'deeptools_input/downreg_H3K9me3.bed'
    script:
        'scripts/parse_DE.py'