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
        down = 'regions/de_down.tsv',
        up = 'regions/de_up.tsv',
        nonde = 'regions/nonde.tsv'
    params:
        padj = config['padj'],
        l2fc = config['l2fc']
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

rule call_H3K4me1_peaks:
    input:
        bam = 'deeptools_input/{sample}.bam',
        ctrl = lambda wildcards: f'deeptools_input/{H3K4me1_CHIPS[wildcards.sample]}.bam'
    output:
        peak = temp('regions/{sample}_peaks.narrowPeak'),
        xls = temp('regions/{sample}_peaks.xls'),
        summits = temp('regions/{sample}_summits.bed')
    shell:'''
    macs3 callpeak -q 1e-2 -t {input.bam} -c {input.ctrl} \
      --keep-dup all \
      --outdir regions \
      -n {wildcards.sample} -f BAMPE -g mm
    '''

rule merge_H3K4me1_peaks:
    input:
        peaks = expand('regions/{sample}_peaks.narrowPeak', sample=H3K4me1_CHIPS.keys())
    output:
        bed = 'regions/H3K4me1.bed',
    shell:'''
    cat {input.peaks} | sort -k1,1 -k2,2n | bedtools merge > {output.bed}
    '''

rule annotate_peaks:
    input:
        gtf = 'deeptools_input/mouse.gtf',
        bed = 'regions/{ann_chip}.bed',
    output:
        beda = temp('regions/{ann_chip}_uropa_allhits.bed'),
        txta = temp('regions/{ann_chip}_uropa_allhits.txt'),
        bedf = temp('regions/{ann_chip}_uropa_finalhits.bed'),
        txtf = 'regions/{ann_chip}_uropa_finalhits.txt',
        json = temp('regions/{ann_chip}_uropa.json'),
        pdf = temp('regions/{ann_chip}_uropa_summary.pdf')
    threads: 10
    shell:'''
    uropa -b {input.bed} -g {input.gtf} --summary \
      --feature gene --distance 100000 100000 \
      --internals 1 -p {wildcards.ann_chip}_uropa -o regions \
      --show-attributes gene_id gene_name
    '''