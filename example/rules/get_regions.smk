def return_peakfiles(chip, samples, broadmarks):
    if chip in broadmarks:
        return expand('regions/{sample}_peaks.broadPeak', sample=[s for s in samples if chip in s])
    else:
        return expand('regions/{sample}_peaks.narrowPeak', sample=[s for s in samples if chip in s])


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
        ctrl = lambda wildcards: f'deeptools_input/{sampleconfig['chipdict'][wildcards.sample]}.bam'
    output:
        xls = temp('regions/{sample}_peaks.xls'),
    params:
        broad = lambda wildcards: '--broad' if any(mark in wildcards.sample for mark in BROADMARKS) else ''
    shell:'''
    macs3 callpeak {params.broad} -q 1e-2 -t {input.bam} -c {input.ctrl} \
      --keep-dup all \
      --outdir regions \
      -n {wildcards.sample} -f BAMPE -g mm
    '''

rule merge_peaks:
    input:
        expand('regions/{sample}_peaks.xls', sample=SAMPLES)
    output:
        bed = 'regions/{chip}.bed',
    params:
        peaks = lambda wildcards: return_peakfiles(wildcards.chip, SAMPLES, BROADMARKS)
    shell:'''
    cat {params.peaks} | sort -k1,1 -k2,2n | bedtools merge > {output.bed}
    '''

rule cleanup_peakfiles:
    input:
        expand('regions/{chip}.bed', chip=CHIPS)
    output:
        temp(touch('regions/peaks_cleaned.txt'))
    shell:'''
    rm -f regions/*_summits.bed
    rm -f regions/*.gappedPeak
    rm -f regions/*.broadPeak
    rm -f regions/*.narrowPeak
    '''

rule annotate_peaks:
    input:
        gtf = 'deeptools_input/mouse.gtf',
        bed = 'regions/{chip}.bed',
        clean = 'regions/peaks_cleaned.txt'
    output:
        beda = temp('regions/{chip}_uropa_allhits.bed'),
        txta = temp('regions/{chip}_uropa_allhits.txt'),
        bedf = temp('regions/{chip}_uropa_finalhits.bed'),
        txtf = 'regions/{chip}_uropa_finalhits.txt',
        json = temp('regions/{chip}_uropa.json'),
        pdf = temp('regions/{chip}_uropa_summary.pdf')
    threads: 10
    shell:'''
    uropa -b {input.bed} -g {input.gtf} --summary \
      --feature gene --distance 100000 100000 \
      --internals 1 -p {wildcards.chip}_uropa -o regions \
      --show-attributes gene_id gene_name
    '''