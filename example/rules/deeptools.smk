rule multibamsummary:
    input:
        bed = 'regions/{mergedpeak}.bed',
        bamfiles = lambda wildcards: expand('deeptools_input/{sample}.bam', sample=[i for i in SAMPLES if wildcards.mergedpeak in i])
    output:
        npz = 'regions/{mergedpeak}_mbs.npz',
        counts = 'regions/{mergedpeak}_mbs.counts'
    threads: 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:'''
    multiBamSummary BED-file -p {threads} -o {output.npz} \
      --BED {input.bed} -b {input.bamfiles} --outRawCounts {output.counts}
    '''

rule de_to_region:
    input:
        down = 'regions/de_down.tsv',
        up = 'regions/de_up.tsv',
        nonde = 'regions/nonde.tsv',
        gtf = 'deeptools_input/mouse.gtf',
    output:
        downbed = 'deeptools_input/downreg_tss.bed',
        upbed = 'deeptools_input/upreg_tss.bed',
        nonbed = 'deeptools_input/nonreg_tss.bed',
        downgtf = 'deeptools_input/downreg_genes.gtf',
        upgtf = 'deeptools_input/upreg_genes.gtf',
        nongtf = 'deeptools_input/nonreg_genes.gtf',
    params:
        l2fc = config['l2fc']
    resources:
        mem_mb = 4000,
        runtime = 1440
    script:
        'scripts/de_to_region.py'


rule parse_de:
    input:
        down = 'regions/de_down.tsv',
        up = 'regions/de_up.tsv',
        nonde = 'regions/nonde.tsv',
        mbs = 'regions/{mergedpeak}_mbs.npz',
        urop = 'regions/{mergedpeak}_uropa_finalhits.txt',
    output:
        down = 'deeptools_input/downreg_{mergedpeak}.bed',
        up = 'deeptools_input/upreg_{mergedpeak}.bed',
        nonde = 'deeptools_input/nonreg_{mergedpeak}.bed'
    params:
        chip = lambda wildcards: wildcards.mergedpeak,
        l2fc = config['l2fc']
    resources:
        mem_mb = 4000,
        runtime = 1440
    script:
        'scripts/parse_DE.py'


rule bamcompare_chip:
    input:
        bam = 'deeptools_input/{chipsample}.bam'
    output:
        bw = 'deeptools_output/chip/{chipsample}.bw'
    threads: 10
    params:
        input = lambda wildcards: f"deeptools_input/{sampleconfig['chipdict'][wildcards.chipsample]}.bam",
        operation = lambda wildcards: 'subtract' if any(mark in wildcards.chipsample for mark in BROADMARKS) else 'log2',
        rar = config['rar'],
        binsize = 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:'''
    bamCompare --operation {params.operation} -b1 {input.bam} -b2 {params.input} \
      -bs {params.binsize} --smoothLength 30 -p {threads} -o {output.bw} --extendReads --centerReads \
      --blackListFileName {params.rar}
    '''

rule bamCoverage_atac:
    input:
        bam = 'deeptools_input/{atacsample}.bam'
    output:
        bw = 'deeptools_output/atac/{atacsample}.bw',
        bw_raw = temp('deeptools_output/atac/{atacsample}_raw.bw'),
        wig = temp('deeptools_output/atac/{atacsample}.wig'),
    params:
        chromsizes = config['chromsizes'],
        rar = config['rar']
    threads: 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:'''
    bamCoverage -b {input.bam} -o {output.bw_raw} \
      --normalizeUsing RPKM -bs 10 --smoothLength 30 -p {threads} \
      --blackListFileName {params.rar}
    wiggletools log 10 {output.bw_raw} > {output.wig}
    wigToBigWig {output.wig} {params.chromsizes} {output.bw}
    '''

rule computeMatrix_chip:
    input:
        bw = expand('deeptools_output/chip/{chipsample}.bw', chipsample=sampleconfig['chipdict'].keys()),
        regions = ["deeptools_input/upreg_{chip}.bed", "deeptools_input/downreg_{chip}.bed", "deeptools_input/nonreg_{chip}.bed"]
    output:
        mat = 'deeptools_output/chip_{chip}.npz'
    params:
        bws = lambda wildcards, input: ' '.join([i for i in input.bw if wildcards.chip in i]),
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw if wildcards.chip in i] ),
        mattype = lambda wildcards: 'scale-regions -a 2000 -b 2000 -m 4000' if wildcards.chip in BROADMARKS else 'reference-point -a 3000 -b 3000 --referencePoint center',
        binsize = 10
    threads: 10
    resources:
        mem_mb = 16000,
        runtime = 1440
    shell:'''
    computeMatrix {params.mattype} -p {threads} \
      -S {params.bws} \
      -o {output.mat} \
      -bs {params.binsize} \
      --missingDataAsZero \
      -R {input.regions} \
      --samplesLabel {params.labels}
    '''

rule plotHeatmap_chip:
    input:
        mat = 'deeptools_output/chip_{chip}.npz'
    output:
        png = 'deeptools_output/chip_{chip}.png'
    params:
        cmap = lambda wildcards: cmap[wildcards.chip]
    resources:
        mem_mb = 4000,
        runtime = 1440
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "\\-3kb" --endLabel "\\+3kb" --colorMap {params.cmap} \
      --xAxisLabel "" \
      --legendLocation none \
      --interpolationMethod bilinear \
      --regionsLabel up down non-de
    '''

rule computeMatrix_atac:
    input:
        bw = expand('deeptools_output/atac/{atacsample}.bw', atacsample=ATACSAMPLES),
        regions = ["deeptools_input/upreg_ATAC.bed", "deeptools_input/downreg_ATAC.bed", "deeptools_input/nonreg_ATAC.bed"]
    output:
        mat = 'deeptools_output/atac.npz'
    params:
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw] )
    threads: 10
    resources:
        mem_mb = 16000,
        runtime = 1440
    shell:'''
    computeMatrix reference-point -p {threads} \
      -S {input.bw} \
      -a 3000 -b 3000 \
      -o {output.mat} \
      --referencePoint center \
      -bs 10 \
      --missingDataAsZero \
      -R {input.regions} \
      --samplesLabel {params.labels}
    '''

rule plotHeatmap_atac:
    input:
        mat = 'deeptools_output/atac.npz'
    output:
        png = 'deeptools_output/atac.png'
    resources:
        mem_mb = 4000,
        runtime = 1440
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel r"-3kb" --endLabel r"+3kb" \
      --colorMap Reds \
      --refPointLabel "TSS" \
      --xAxisLabel "" \
      --interpolationMethod bilinear \
      --regionsLabel up down non-de
    '''

rule computeMatrix_meth:
    input:
        bw = expand('deeptools_input/{bssample}_CpG.bw', bssample=BSSAMPLES),
        regions = ["regions/meth_up.bed", "regions/meth_down.bed", "regions/meth_nonde.bed"]
    output:
        mat = 'deeptools_output/meth.npz'
    threads: 10
    resources:
        mem_mb = 16000,
        runtime = 1440
    params:
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw] )
    shell:'''
    computeMatrix reference-point -p {threads} \
      -S {input.bw} \
      -a 500 -b 500 \
      -o {output.mat} \
      --referencePoint center \
      -bs 1 \
      -R {input.regions} \
      --samplesLabel {params.labels}
    '''

rule plotHeatmap_meth:
    input:
        mat = 'deeptools_output/meth.npz'
    output:
        png = 'deeptools_output/meth.png'
    resources:
        mem_mb = 4000,
        runtime = 1440
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "\\-0.5kb" --endLabel "\\+0.5kb" \
      --colorMap Grays_r \
      --interpolationMethod bilinear \
      --regionsLabel up down non-de \
      --legendLocation upper-left
    '''
