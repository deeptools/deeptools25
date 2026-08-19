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

rule bamCoverage_RNA:
    input:
        bam = 'deeptools_input/{rnasample}.bam'
    output:
        bw = 'deeptools_output/rna/{rnasample}.bw',
    params:
        chromsizes = config['chromsizes'],
        rar = config['rar']
    threads: 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:'''
    bamCoverage -b {input.bam} -o {output.bw} \
      --normalizeUsing BPM -bs 10 --smoothLength 30 -p {threads} \
      --blackListFileName {params.rar}
    '''

rule computeMatrix_rna:
    input:
        bw = expand('deeptools_output/rna/{rnasample}.bw', rnasample=RNASAMPLES),
        regions = ["deeptools_input/upreg_genes.gtf", "deeptools_input/downreg_genes.gtf", "deeptools_input/nonreg_genes.gtf"]
    output:
        mat = 'deeptools_output/rna.npz'
    params:
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw] )
    threads: 10
    resources:
        mem_mb = 16000,
        runtime = 1440
    shell:'''
    computeMatrix scale-regions -p {threads} \
      -S {input.bw} \
      -a 200 -b 200 \
      -o {output.mat} \
      -bs 10 \
      --missingDataAsZero \
      -R {input.regions} \
      --samplesLabel {params.labels} \
      --metagene
    '''

rule plotHeatmap_rna:
    input:
        mat = 'deeptools_output/rna.npz'
    output:
        png = 'deeptools_output/rna.png'
    resources:
        mem_mb = 4000,
        runtime = 1440
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "TSS" --endLabel "TES" \
      --colorMap Blues \
      --regionsLabel up down non-de \
      --zMax 2
    '''

def meth_zmax(range_file):
    with open(range_file) as fh:
        rows = [l.rstrip('\n').split('\t') for l in fh if l.strip()]
    header, vals = rows[0], rows[1:]
    idx = next((i for i, h in enumerate(header) if h.strip().lower() == 'max'), None)
    if idx is None:
        idx = next((i for i, h in enumerate(header) if '90' in h), len(header) - 1)
    return round(max(float(row[idx]) for row in vals) * 1.1, 3)

rule meth_bins:
    input:
        down = 'regions/de_down.tsv',
        up = 'regions/de_up.tsv',
        nonde = 'regions/nonde.tsv',
        gtf = 'deeptools_input/mouse.gtf'
    output:
        bed = 'deeptools_input/meth_bins.bed'
    params:
        binsize = 50,
        flank = 3000
    resources:
        mem_mb = 2000,
        runtime = 1440
    script:
        'scripts/make_meth_bins.py'

rule meth_bigwigsummary:
    input:
        bed = 'deeptools_input/meth_bins.bed',
        bw = expand('deeptools_input/{bssample}_CpG.bw', bssample=BSSAMPLES)
    output:
        npz = 'regions/meth_bins.npz',
        counts = 'regions/meth_bins.counts'
    threads: 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:'''
    multiBigwigSummary BED-file -p {threads} -o {output.npz} \
      --BED {input.bed} -b {input.bw} --smartLabels --outRawCounts {output.counts}
    '''

rule parse_meth:
    input:
        down = 'regions/de_down.tsv',
        up = 'regions/de_up.tsv',
        nonde = 'regions/nonde.tsv',
        bins = 'deeptools_input/meth_bins.bed',
        counts = 'regions/meth_bins.counts'
    output:
        down = 'deeptools_input/downreg_meth.bed',
        up = 'deeptools_input/upreg_meth.bed',
        nonde = 'deeptools_input/nonreg_meth.bed'
    params:
        mdiff = config['mdiff']
    resources:
        mem_mb = 4000,
        runtime = 1440
    script:
        'scripts/parse_meth.py'

rule computeMatrix_meth:
    input:
        bw = expand('deeptools_input/{bssample}_CpG.bw', bssample=BSSAMPLES),
        regions = ["deeptools_input/upreg_meth.bed", "deeptools_input/downreg_meth.bed", "deeptools_input/nonreg_meth.bed"]
    output:
        mat = 'deeptools_output/meth.npz'
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
      -bs 50 \
      --missingDataAsZero \
      -R {input.regions} \
      --samplesLabel {params.labels}
    '''

rule meth_datarange:
    input:
        mat = 'deeptools_output/meth.npz'
    output:
        txt = 'deeptools_output/meth_datarange.txt'
    resources:
        mem_mb = 2000,
        runtime = 1440
    shell:'''
    computeMatrixOperations dataRange -m {input.mat} > {output.txt}
    '''

rule plotHeatmap_meth:
    input:
        mat = 'deeptools_output/meth.npz',
        range = 'deeptools_output/meth_datarange.txt'
    output:
        png = 'deeptools_output/meth.png'
    params:
        zmax = lambda wildcards, input: meth_zmax(input.range)
    resources:
        mem_mb = 4000,
        runtime = 1440
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "TSS" --endLabel "TES" --colorMap YlOrRd --zMin 0 --zMax {params.zmax} \
      --interpolationMethod bilinear \
      --regionsLabel up down non-de \
      --sortRegions descend
    '''
