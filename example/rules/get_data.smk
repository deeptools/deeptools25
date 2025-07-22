
if config['source'] == 'raw':
    ITS = [(k,v,'fq') for k,v in sampleconfig['samples'].items()]

    rule fqfiles:
        output:
            expand('fq/{sample}_{R}.fastq.gz', sample=SAMPLES, R=['R1', 'R2'])
        params:
            its = ITS
        threads: 10
        script:
            'scripts/get_data.py'

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

elif config['source'] == 'zenodo':
    print("Zenodo mode")





# Zenodo mode
