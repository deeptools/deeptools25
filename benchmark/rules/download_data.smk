rule download_all:
    output:
        expand("zenodo/{f}", f=ALLFILES)
    params:
        odir="zenodo",
        zenodo_id=sampleconfig["zenodo"]["ID"],
    threads: 10
    script:
        "scripts/download_zenodo.py"


rule validate_files:
    input:
        "zenodo/{file}"
    output:
        touch("zenodo/{file}.valid")
    params:
        exp=lambda wildcards: sampleconfig["md5sums"][wildcards.file],
    run:
        import hashlib
        hash_md5 = hashlib.md5()
        with open(input[0], "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        got = hash_md5.hexdigest()
        assert got == params.exp, f"MD5 mismatch for {input[0]}: expected {params.exp}, got {got}"


rule cram_to_bam:
    input:
        cram=lambda wc: "zenodo/{cramfile}.cram".format(cramfile=wc.cramfile),
        valid="zenodo/{cramfile}.cram.valid",
        fna=lambda wc: f"zenodo/{'human' if 'human' in wc.cramfile else 'triticum'}.fna",
        fna_valid=lambda wc: f"zenodo/{'human' if 'human' in wc.cramfile else 'triticum'}.fna.valid",
    output:
        bam="bamfiles/{cramfile}.bam",
        bai="bamfiles/{cramfile}.bam.bai",
    threads: 10
    shell:
        "samtools view -@ {threads} -T {input.fna} -b -o {output.bam} {input.cram} && "
        "samtools index -@ {threads} {output.bam}"
