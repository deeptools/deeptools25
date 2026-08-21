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

wildcard_constraints:
    condition = 'ctrl|ko'

rule bigwigAverage_meth:
    input:
        bw = lambda wildcards: expand(
            'deeptools_input/{bssample}_CpG.bw',
            bssample=[s for s in BSSAMPLES if f'_{wildcards.condition}_' in s]
        )
    output:
        bw = 'deeptools_output/meth_{condition}_avg.bw'
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:'''
    bigwigAverage -b {input.bw} -p {threads} -o {output.bw}
    '''

rule meth_diff:
    input:
        ko = 'deeptools_output/meth_ko_avg.bw',
        ctrl = 'deeptools_output/meth_ctrl_avg.bw'
    output:
        bw = 'deeptools_output/meth_diff.bw'
    resources:
        mem_mb = 4000,
        runtime = 1440
    shell:'''
    bigwigCompare -b1 {input.ko} -b2 {input.ctrl} --operation subtract -bs 50 -o {output.bw}
    '''

rule computeMatrix_meth_diff:
    input:
        bw = 'deeptools_output/meth_diff.bw',
        regions = ["deeptools_input/upreg_meth.bed", "deeptools_input/downreg_meth.bed", "deeptools_input/nonreg_meth.bed"]
    output:
        mat = 'deeptools_output/meth_diff.npz'
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
      -R {input.regions} \
      --samplesLabel "ko - ctrl"
    '''

rule meth_diff_datarange:
    input:
        mat = 'deeptools_output/meth_diff.npz'
    output:
        txt = 'deeptools_output/meth_diff_datarange.txt'
    resources:
        mem_mb = 2000,
        runtime = 1440
    shell:'''
    computeMatrixOperations dataRange -m {input.mat} > {output.txt}
    '''

def meth_diff_zmax(range_file):
    # zMax is bounded by the 10th/90th percentile (not the raw max) so a
    # handful of outlier bins can't blow out the whole color scale.
    with open(range_file) as fh:
        rows = [l.rstrip('\n').split('\t') for l in fh if l.strip()]
    header, vals = rows[0], rows[1:]
    idx10 = next(i for i, h in enumerate(header) if '10' in h)
    idx90 = next(i for i, h in enumerate(header) if '90' in h)
    bound = max(max(abs(float(row[idx10])), abs(float(row[idx90]))) for row in vals)
    return round(max(bound, 1e-6) * 1.1, 3)

rule plotHeatmap_meth_diff:
    input:
        mat = 'deeptools_output/meth_diff.npz',
        range = 'deeptools_output/meth_diff_datarange.txt'
    output:
        png = 'deeptools_output/meth_diff.png'
    params:
        zmax = lambda wildcards, input: meth_diff_zmax(input.range)
    resources:
        mem_mb = 4000,
        runtime = 1440
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "\\-3kb" --endLabel "\\+3kb" \
      --colorMap RdBu_r --zMin=-{params.zmax} --zMax={params.zmax} \
      --interpolationMethod bilinear \
      --regionsLabel up down non-de \
      --sortUsing max \
      --legendLocation upper-left
    '''
