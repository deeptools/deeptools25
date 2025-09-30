def return_peakfiles(mergedpeak, samples, broadmarks):
    if mergedpeak in broadmarks:
        return expand('regions/{sample}_peaks.broadPeak', sample=[s for s in samples if mergedpeak in s])
    else:
        return expand('regions/{sample}_peaks.narrowPeak', sample=[s for s in samples if mergedpeak in s])


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
        bam = expand('deeptools_input/{sample}.bam', sample=SAMPLES),
    output:
        xls = temp('regions/{peaksample}_peaks.xls'),
    params:
        bam = lambda wildcards: f'deeptools_input/{wildcards.peaksample}.bam',
        ctrl = lambda wildcards: f'-c deeptools_input/{sampleconfig['chipdict'][wildcards.peaksample]}.bam' if wildcards.peaksample in sampleconfig['chipdict'] else '',
        broad = lambda wildcards: '--broad' if any(mark in wildcards.peaksample for mark in BROADMARKS) else '',
        atacpar = lambda wildcards: '--nomodel --shift -75 --extsize 150' if 'ATAC' in wildcards.peaksample else '',
    shell:'''
    macs3 callpeak {params.broad} -q 1e-2 -t {params.bam} {params.ctrl} \
      --keep-dup all \
      --outdir regions \
      {params.atacpar} -n {wildcards.peaksample} -f BAMPE -g mm
    '''

rule merge_peaks:
    input:
        expand('regions/{peaksample}_peaks.xls', peaksample = ATACSAMPLES + list(sampleconfig['chipdict'].keys()))
    output:
        bed = 'regions/{mergedpeak}.bed',
    params:
        peaks = lambda wildcards: return_peakfiles(wildcards.mergedpeak, SAMPLES, BROADMARKS)
    shell:'''
    cat {params.peaks} | sort -k1,1 -k2,2n | bedtools merge > {output.bed}
    '''

rule cleanup_peakfiles:
    input:
        expand('regions/{mergedpeak}.bed', mergedpeak = ['ATAC'] + CHIPS)
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
        bed = 'regions/{mergedpeak}.bed',
        clean = 'regions/peaks_cleaned.txt'
    output:
        beda = temp('regions/{mergedpeak}_uropa_allhits.bed'),
        txta = temp('regions/{mergedpeak}_uropa_allhits.txt'),
        bedf = temp('regions/{mergedpeak}_uropa_finalhits.bed'),
        txtf = 'regions/{mergedpeak}_uropa_finalhits.txt',
        json = temp('regions/{mergedpeak}_uropa.json'),
        pdf = temp('regions/{mergedpeak}_uropa_summary.pdf')
    threads: 10
    shell:'''
    uropa -b {input.bed} -g {input.gtf} --summary \
      --feature gene --distance 100000 100000 \
      --internals 1 -p {wildcards.mergedpeak}_uropa -o regions \
      --show-attributes gene_id gene_name
    '''