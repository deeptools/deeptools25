rule bamcompare_chip:
    input:
        bam = 'deeptools_input/{chipsample}.bam'
    output:
        bw = 'deeptools_output/chip/{chipsample}.bw'
    threads: 10
    params:
        input = lambda wildcards: f'deeptools_input/{sampleconfig['chipdict'][wildcards.chipsample]}.bam'
    shell:'''
    bamCompare -b1 {input.bam} -b2 {params.input} \
      -bs 10 -p {threads} -o {output.bw}
    '''

rule bamCoverage_atac:
    input:
        bam = 'deeptools_input/{atacsample}.bam'
    output:
        bw = 'deeptools_output/atac/{atacsample}.bw'
    threads: 10
    shell:'''
    bamCoverage -b {input.bam} -o {output.bw} \
      --normalizeUsing RPKM -bs 10 -p {threads}
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
      --normalizeUsing RPKM -bs 10 -p {threads}
    wiggletools log 10 {output.bw_raw} > {output.wig}
    wigToBigWig {output.wig} {params.chromsizes} {output.bw}
    '''

rule computeMatrix_chip:
    input:
        bw = expand('deeptools_output/chip/{chipsample}.bw', chipsample=sampleconfig['chipdict'].keys())
    output:
        mat = 'deeptools_output/chip_{chip}.npz'
    params:
        bws = lambda wildcards, input: ' '.join([i for i in input.bw if wildcards.chip in i]),
        labels = lambda wildcards, input: ' '.join( [i.split('_')[3] + '-' + i.split('_')[4] for i in input.bw if wildcards.chip in i] )
    threads: 10
    shell:'''
    computeMatrix reference-point -p {threads} \
      -S {params.bws} \
      -a 3000 -b 3000 \
      -o {output.mat} \
      --referencePoint center \
      -bs 10 \
      --missingDataAsZero \
      -R /data/manke/processing/deboutte/tmp/region_generator/uptss.bed /data/manke/processing/deboutte/tmp/region_generator/downtss.bed \
      --samplesLabel {params.labels}
    '''

rule plotHeatmap_chip:
    input:
        mat = 'deeptools_output/chip_{chip}.npz'
    output:
        png = 'deeptools_output/chip_{chip}.png'
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel "\\-3kb" --endLabel "\\+3kb" --interpolationMethod bicubic
    '''

rule computeMatrix_atac:
    input:
        bw = expand('deeptools_output/atac/{atacsample}.bw', atacsample=ATACSAMPLES)
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
      -bs 10 \
      --missingDataAsZero \
      -R /data/manke/processing/deboutte/tmp/region_generator/uptss.bed /data/manke/processing/deboutte/tmp/region_generator/downtss.bed \
      --samplesLabel {params.labels}
    '''

rule plotHeatmap_atac:
    input:
        mat = 'deeptools_output/atac.npz'
    output:
        png = 'deeptools_output/atac.png'
    shell:'''
    plotHeatmap -m {input.mat} -out {output.png} \
      --startLabel r"-3kb" --endLabel r"+3kb" --interpolationMethod bicubic
    '''

rule computeMatrix_rna:
    input:
        bw = expand('deeptools_output/rna/{rnasample}.bw', rnasample=RNASAMPLES)
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
      -bs 10 \
      --missingDataAsZero \
      -R /data/manke/processing/deboutte/tmp/region_generator/upreg.gtf /data/manke/processing/deboutte/tmp/region_generator/downreg.gtf \
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
      --startLabel "TSS" --endLabel "TES" --interpolationMethod bicubic
    '''