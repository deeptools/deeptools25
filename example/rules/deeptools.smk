rule multibamsummary:
    input:
        bed = 'regions/{ann_chip}.bed',
        bamfiles = lambda wildcards: expand('deeptools_input/{sample}.bam', sample=[i for i in SAMPLES if wildcards.ann_chip in i])
    output:
        npz = 'regions/{ann_chip}_mbs.npz'
    threads: 10
    shell:'''
    multiBamSummary BED-file -p {threads} -o {output.npz} \
      --BED {input.bed} -b {input.bamfiles}
    '''
rule parse_regions:
    input:
        down = 'regions/de_down.tsv',
        up = 'regions/de_up.tsv',
        nonde = 'regions/nonde.tsv',
        gtf = 'deeptools_input/mouse.gtf',
        uro_k27 = 'regions/H3K27me3_uropa_finalhits.txt',
        uro_k9 = 'regions/H3K9me3_uropa_finalhits.txt',
        uro_k4me1 = 'regions/H3K4me1_uropa_finalhits.txt',
        mbs_k27 = 'regions/H3K27me3_mbs.npz',
        mbs_k9 = 'regions/H3K9me3_mbs.npz',
        mbs_k4me1 = 'regions/H3K4me1_mbs.npz'
    output:
        downbed = 'deeptools_input/downreg_tss.bed',
        upbed = 'deeptools_input/upreg_tss.bed',
        nonbed = 'deeptools_input/nonreg_tss.bed',
        downgtf = 'deeptools_input/downreg_genes.gtf',
        upgtf = 'deeptools_input/upreg_genes.gtf',
        nongtf = 'deeptools_input/nonreg_genes.gtf',
        k27_down = 'deeptools_input/downreg_H3K27me3.bed',
        k9_down = 'deeptools_input/downreg_H3K9me3.bed',
        k4_down = 'deeptools_input/downreg_H3K4me1.bed'
    params:
        l2fc = config['l2fc']
    script:
        'scripts/parse_DE.py'

rule bamcompare_chip:
    input:
        bam = 'deeptools_input/{chipsample}.bam'
    output:
        bw = 'deeptools_output/chip/{chipsample}.bw'
    threads: 10
    params:
        input = lambda wildcards: f'deeptools_input/{sampleconfig['chipdict'][wildcards.chipsample]}.bam',
        bs = lambda wildcards: '-bs 100' if 'H3K27me3' in wildcards.chipsample or 'H3K9me3' in wildcards.chipsample else '-bs 10',
    shell:'''
    bamCompare -b1 {input.bam} -b2 {params.input} \
      {params.bs} -p {threads} -o {output.bw} --extendReads
    '''

rule bamCoverage_atac:
    input:
        bam = 'deeptools_input/{atacsample}.bam'
    output:
        bw = 'deeptools_output/atac/{atacsample}.bw',
    params:
        chromsizes = config['chromsizes']
    threads: 10
    shell:'''
    bamCoverage -b {input.bam} -o {output.bw} \
      --normalizeUsing RPKM -bs 50 -p {threads}
    '''

rule bamCoverage_RNA:
    input:
        bam = 'deeptools_input/{rnasample}.bam'
    output:
        bw = 'deeptools_output/rna/{rnasample}.bw',
        bw_raw = temp('deeptools_output/rna/{rnasample}_raw.bw'),
        wig = temp('deeptools_output/rna/{rnasample}.wig'),
    params:
        chromsizes = config['chromsizes']
    threads: 10
    shell:'''
    bamCoverage -b {input.bam} -o {output.bw_raw} \
      --normalizeUsing RPKM -bs 50 -p {threads}
    wiggletools log 10 {output.bw_raw} > {output.wig}
    wigToBigWig {output.wig} {params.chromsizes} {output.bw}
    '''

rule computeMatrix_chip:
    input:
        bw = expand('deeptools_output/chip/{chipsample}.bw', chipsample=sampleconfig['chipdict'].keys()),
        regions = lambda wildcards: (
            ["deeptools_input/downreg_tss.bed", "deeptools_input/upreg_tss.bed", "deeptools_input/nonreg_tss.bed"] if wildcards.chip in ['H3K4me3','H3K27ac']
            else ["deeptools_input/downreg_H3K4me1.bed"] if wildcards.chip == 'H3K4me1'
            else ["deeptools_input/downreg_H3K27me3.bed"] if wildcards.chip == 'H3K27me3'
            else ["deeptools_input/downreg_H3K9me3.bed"])
    output:
        mat = 'deeptools_output/chip_{chip}.npz'
    params:
        bws = lambda wildcards, input: ' '.join([i for i in input.bw if wildcards.chip in i]),
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw if wildcards.chip in i] ),
        mattype = lambda wildcards: 'scale-regions -a 200 -b 200' if wildcards.chip in INH_CHIP else 'reference-point -a 3000 -b 3000 --referencePoint center',
        bs = lambda wildcards: '-bs 100' if 'H3K27me3' == wildcards.chip or 'H3K9me3' == wildcards.chip else '-bs 10',
    threads: 10
    shell:'''
    computeMatrix {params.mattype} -p {threads} \
      -S {params.bws} \
      -o {output.mat} \
      {params.bs} \
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
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "\\-3kb" --endLabel "\\+3kb" --colorMap {params.cmap} \
      --xAxisLabel "" \
      --legendLocation none
    '''

rule computeMatrix_atac:
    input:
        bw = expand('deeptools_output/atac/{atacsample}.bw', atacsample=ATACSAMPLES),
        regions = ["deeptools_input/downreg_tss.bed", "deeptools_input/upreg_tss.bed", "deeptools_input/nonreg_tss.bed"]
    output:
        mat = 'deeptools_output/atac.npz'
    params:
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw] )
    threads: 10
    shell:'''
    computeMatrix reference-point -p {threads} \
      -S {input.bw} \
      -a 3000 -b 3000 \
      -o {output.mat} \
      --referencePoint center \
      -bs 50 \
      --missingDataAsZero \
      -R {input.regions} \
      --samplesLabel {params.labels} \
      --scale 2
    '''

rule plotHeatmap_atac:
    input:
        mat = 'deeptools_output/atac.npz'
    output:
        png = 'deeptools_output/atac.png'
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel r"-3kb" --endLabel r"+3kb" --colorMap Blues \
      --refPointLabel "TSS" \
      --regionsLabel "downreg" "upreg" "nonde" \
      --xAxisLabel "" \
      --legendLocation none
    '''

rule computeMatrix_rna:
    input:
        bw = expand('deeptools_output/rna/{rnasample}.bw', rnasample=RNASAMPLES),
        regions = ["deeptools_input/downreg_genes.gtf", "deeptools_input/upreg_genes.gtf", "deeptools_input/nonreg_genes.gtf"]
    output:
        mat = 'deeptools_output/rna.npz'
    params:
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw] )
    threads: 10
    shell:'''
    computeMatrix scale-regions -p {threads} \
      -S {input.bw} \
      -a 200 -b 200 \
      -o {output.mat} \
      -bs 50 \
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
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "TSS" --endLabel "TES" --colorMap Reds
    '''

rule computeMatrix_meth:
    input:
        bw = expand('deeptools_input/{bssample}_CpG.bw', bssample=BSSAMPLES),
        regions = ["deeptools_input/downreg_tss.bed", "deeptools_input/upreg_tss.bed", "deeptools_input/nonreg_tss.bed"]
    output:
        mat = 'deeptools_output/meth.npz'
    params:
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw] )
    threads: 10
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

rule plotHeatmap_meth:
    input:
        mat = 'deeptools_output/meth.npz'
    output:
        png = 'deeptools_output/meth.png'
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "TSS" --endLabel "TES" --colorMap Greys --zMin 0
    '''