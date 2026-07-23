rule download_cram:
    output:
        expand('zenodo_dl/{sample}.cram', sample=SAMPLES),
        fna = 'zenodo_dl/mouse.fna',
        gtf = 'zenodo_dl/mouse.gtf'
    params:
        odir = 'zenodo_dl',
        zenodo_id = sampleconfig['zenodo']['ID'],
        only_fna = False
    threads: 10
    script:
        'scripts/download_zenodo.py'

rule validate_cramfiles:
    input:
        cramfile = 'zenodo_dl/{sample}.cram',
    output:
        touch('zenodo_dl/{sample}.valid')
    params:
        exp = lambda wildcards: sampleconfig['md5sums']['cram'][wildcards.sample],
    run:
        import hashlib
        hash_md5 = hashlib.md5()
        with open(input.cramfile, 'rb') as f:
            for _ in iter(lambda: f.read(4096), b""):
                hash_md5.update(_)
        assert hash_md5.hexdigest() == params.exp, f"MD5 mismatch for {input.cramfile}, expected {params.exp}, got {hash_md5.hexdigest()}"

rule prep_deeptools_input:
    input:
        cramfile = 'zenodo_dl/{sample}.cram',
        valid = 'zenodo_dl/{sample}.valid',
        fna = 'zenodo_dl/mouse.fna'
    output:
        bam = 'deeptools_input/{sample}.bam',
        bai = 'deeptools_input/{sample}.bam.bai'
    threads: 10
    run:
        shell('samtools view -f 0x2 -@ {threads} -T {input.fna} -b -o {output.bam} {input.cramfile}')
        shell('samtools index -@ {threads} {output.bam}')

rule generate_bs_bedgraph_zenodo:
    input:
        bam = 'deeptools_input/{bssample}.bam',
        bai = 'deeptools_input/{bssample}.bam.bai',
        fna = 'zenodo_dl/mouse.fna'
    output:
        bg = temp('deeptools_input/{bssample}_CpG.bedGraph'),
        bgs = temp('deeptools_input/{bssample}_CpG_subset.bedGraph'),
        bw = 'deeptools_input/{bssample}_CpG.bw'
    params:
        chromsizes = config['chromsizes']
    threads: 10
    shell:'''
    MethylDackel extract -@ {threads} {input.fna} {input.bam}
    cut -f1,2,3,4 {output.bg} > {output.bgs}
    bedGraphToBigWig {output.bgs} {params.chromsizes} {output.bw}
    '''

rule ship_zen_fna_gtf:
    input:
        fna = 'zenodo_dl/mouse.fna',
        gtf = 'zenodo_dl/mouse.gtf'
    output:
        fna = 'deeptools_input/mouse.fna',
        gtf = 'deeptools_input/mouse.gtf'
    run:
        import shutil
        shutil.copy2(input.fna, output.fna)
        shutil.copy2(input.gtf, output.gtf)
