rule download_data:
    output:
        ofile = "zenodo/{file}",
        valid = "zenodo/{file}.valid"
    params:
        odir="zenodo",
        zenodo_id=sampleconfig["zenodo"]["ID"],
    threads: 2
    resources:
      mem_mb = 10000,
      runtime = 1440
    script:
        "scripts/download_zenodo.py"

rule cram_to_bam:
    input:
        cram=lambda wc: "zenodo/{cramfile}.cram".format(cramfile=wc.cramfile),
        fna=lambda wc: f"zenodo/{'human' if 'human' in wc.cramfile else 'triticum'}.fna",
    output:
        bam="bamfiles/{cramfile}.bam"
    params:
        ix_param = lambda wc: '-c' if 'triticum' in wc.cramfile else ''
    threads: 10
    resources:
      mem_mb = 10000,
      runtime = 1440
    shell:"""
    samtools view -@ {threads} -T {input.fna} -b -o {output.bam} {input.cram}
    samtools index -@ {threads} {params.ix_param} {output.bam}
    """
