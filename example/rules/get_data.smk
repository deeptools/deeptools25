
if config['source'] == 'raw':
    ITS = [(k,v,'fq') for k,v in sampleconfig['samples'].items()]

    rule fqfiles:
        output:
            expand('fq/{sample}_{R}.fastq.gz', sample=SAMPLES, R=['R1', 'R2'])
        params:
            its = ITS
        threads: 10
        script:
            'scripts/download_fq.py'
    
    rule download_fna:
        output:
            fna = 'fq/mouse.fna'
        params:
            odir = 'fq',
            zenodo_id = sampleconfig['zenodo']['ID'],
            only_fna = True
        threads: 10
        script:
            'scripts/download_zenodo.py'

    rule validate_fqfiles:
        input:
            r1 = 'fq/{sample}_R1.fastq.gz',
            r2 = 'fq/{sample}_R2.fastq.gz'
        output:
            touch('fq/{sample}.valid')
        params:
            exp_r1 = lambda wildcards: sampleconfig['md5sums']['fq'][wildcards.sample][0],
            exp_r2 = lambda wildcards: sampleconfig['md5sums']['fq'][wildcards.sample][1],
        run:
            import hashlib
            for input_file, exp in zip([input.r1, input.r2], [params.exp_r1, params.exp_r2]):
                hash_md5 = hashlib.md5()
                with open(input_file, 'rb') as f:
                    for _ in iter(lambda: f.read(4096), b""):
                        hash_md5.update(_)
                assert hash_md5.hexdigest() == exp, f"MD5 mismatch for {input_file}, expected {exp}, got {hash_md5.hexdigest()}"

    rule prep_snakepipes:
        input:
            expand('fq/{sample}.valid', sample=SAMPLES)
        output:
            atacd = directory('fq/atac'),
            bsd = directory('fq/bs'),
            rnad = directory('fq/rna'),
            chipd = directory('fq/chip')
        run:
            import os
            import glob
            os.makedirs('fq/atac', exist_ok=True)
            os.makedirs('fq/bs', exist_ok=True)
            os.makedirs('fq/rna', exist_ok=True)
            os.makedirs('fq/chip', exist_ok=True)
            # Link fq files to dirs
            ## ATAC
            for _ in glob.glob(f"fq/*ATAC*fastq.gz"):
                os.symlink('../' + _.replace('fq/', ''), f"fq/atac/{os.path.basename(_)}")
            ## BS
            for _ in glob.glob(f"fq/*BS*fastq.gz"):
                os.symlink('../' + _.replace('fq/', ''), f"fq/bs/{os.path.basename(_)}")
            ## RNA
            for _ in glob.glob(f"fq/*RNA*fastq.gz"):
                os.symlink('../' + _.replace('fq/', ''), f"fq/rna/{os.path.basename(_)}")
            ## ChIP
            for _ in glob.glob(f"fq/*.fastq.gz"):
                if 'ATAC' not in _ and 'BS' not in _ and 'RNA' not in _:
                    os.symlink('../' + _.replace('fq/', ''), f"fq/chip/{os.path.basename(_)}")
    
    rule run_snakepipes_atac:
        localrule: True
        input:
            atacd = 'fq/atac'
        output:
            of = 'snakePipes/ATAC/DNAmapping_snakePipes.done',
            od = directory('snakePipes/ATAC/')
        conda:
            config['snakepipes_env']
        shell:'''
        DNAmapping --dedup --mapq 3 --fastqc -i {input.atacd} -o snakePipes/ATAC mm39_ens106
        '''
    
    rule run_snakepipes_rna:
        localrule: True
        input:
            rnad = 'fq/rna',
            atac = 'snakePipes/ATAC/DNAmapping_snakePipes.done'
        output:
            of = 'snakePipes/RNA/mRNAseq_snakePipes.done',
            od = directory('snakePipes/RNA/')
        conda:
            config['snakepipes_env']
        shell:'''
        mRNAseq -i {input.rnad} -o snakePipes/RNA mm39_ens106
        '''
    
    rule run_snakepipes_bs:
        localrule: True
        input:
            bsd = 'fq/bs',
            atac = 'snakePipes/ATAC/DNAmapping_snakePipes.done',
            rna = 'snakePipes/RNA/mRNAseq_snakePipes.done',
        output:
            of = 'snakePipes/BS/WGBS_snakePipes.done',
            od = directory('snakePipes/BS/')
        conda:
            config['snakepipes_env']
        shell:'''
        WGBS -i {input.bsd} -o snakePipes/BS mm39_ens106
        '''
    
    rule run_snakepipes_chip:
        localrule: True
        input:
            chipd = 'fq/chip',
            atac = 'snakePipes/ATAC/DNAmapping_snakePipes.done',
            rna = 'snakePipes/RNA/mRNAseq_snakePipes.done',
            bs = 'snakePipes/BS/WGBS_snakePipes.done',
        output:
            of = 'snakePipes/ChIP/DNAmapping_snakePipes.done',
            od = directory('snakePipes/ChIP/')
        conda:
            config['snakepipes_env']
        shell:'''
        DNAmapping --dedup --mapq 3 --fastqc  -i {input.chipd} -o snakePipes/ChIP mm39_ens106
        '''
    
    rule prep_deeptools_input:
        input:
            atacd = 'snakePipes/ATAC',
            rnad = 'snakePipes/RNA',
            bsd = 'snakePipes/BS',
            chipd = 'snakePipes/ChIP'
        output:
            expand("deeptools_input/{sample}.bam", sample=SAMPLES)
        run:
            import glob
            import shutil
            from pathlib import Path
            # atac
            for bdir in [input.atacd, input.rnad, input.bsd, input.chipd]:
                for bam in (Path(bdir) / 'filtered_bam').glob('*.bam'):
                    _of = Path('deeptools_input') / bam.name.replace('.filtered.bam', '.bam')
                    shutil.copy2(bam, _of)
                    # ship bai too.
                    bai = bam.parent / (bam.name + '.bai')
                    _ofbai = Path('deeptools_input') / bai.name.replace('.filtered.bam', '.bam')
                    shutil.copy2(bai, _ofbai)
    
    rule generate_bs_bedgraph:
        input:
            bam = 'deeptools_input/{bssample}.bam',
            fna = 'fq/mouse.fna'
        output:
            'deeptools_input/{bssample}_CpG.bedGraph'
        threads: 10
        shell:'''
        MethylDackel extract -@ {threads} {input.fna} {input.bam}
        '''

elif config['source'] == 'zenodo':
    rule download_cram:
        output:
            expand('zenodo_dl/{sample}.cram', sample=SAMPLES),
            fna = 'zenodo_dl/mouse.fna'
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
            'deeptools_input/{bssample}_CpG.bedGraph'
        threads: 10
        shell:'''
        MethylDackel extract -@ {threads} {input.fna} {input.bam}
        '''