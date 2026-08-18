rule download_data:
    output:
        ofile = "zenodo/{file}",
        valid = "zenodo/{file}.valid"
    params:
        odir="zenodo",
        zenodo_id=sampleconfig["zenodo"]["ID"],
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 1440
    script:
        "scripts/download_zenodo.py"


rule prep_deeptools_input:
    input:
        cramfile = 'zenodo/{sample}.cram',
        valid = 'zenodo/{sample}.cram.valid',
        fna = 'zenodo/mouse.fna'
    output:
        bam = 'deeptools_input/{sample}.bam',
        bai = 'deeptools_input/{sample}.bam.bai'
    threads: 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    run:
        shell('samtools view -f 0x2 -@ {threads} -T {input.fna} -b -o {output.bam} {input.cramfile}')
        shell('samtools index -@ {threads} {output.bam}')

rule generate_bs_bedgraph_zenodo:
    input:
        bam = 'deeptools_input/{bssample}.bam',
        bai = 'deeptools_input/{bssample}.bam.bai',
        fna = 'zenodo/mouse.fna'
    output:
        bg = temp('deeptools_input/{bssample}_CpG.bedGraph'),
        bgs = temp('deeptools_input/{bssample}_CpG_subset.bedGraph'),
        bw = 'deeptools_input/{bssample}_CpG.bw'
    params:
        chromsizes = config['chromsizes']
    threads: 10
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:'''
    MethylDackel extract -@ {threads} {input.fna} {input.bam}
    cut -f1,2,3,4 {output.bg} > {output.bgs}
    bedGraphToBigWig {output.bgs} {params.chromsizes} {output.bw}
    '''

rule ship_zen_fna_gtf:
    input:
        fna = 'zenodo/mouse.fna',
        gtf = 'zenodo/mouse.gtf'
    output:
        fna = 'deeptools_input/mouse.fna',
        gtf = 'deeptools_input/mouse.gtf'
    resources:
        mem_mb = 2000,
        runtime = 1440
    run:
        import shutil
        shutil.copy2(input.fna, output.fna)
        shutil.copy2(input.gtf, output.gtf)
