rule download_data:
    output:
        ofile = "zenodo/{file}",
        valid = "zenodo/{file}.valid"
    params:
        odir="zenodo",
        zenodo_id=sampleconfig["zenodo"]["ID"],
    threads: 2
    script:
        "scripts/download_zenodo.py"

rule cram_to_bam:
    input:
        cram=lambda wc: "zenodo/{cramfile}.cram".format(cramfile=wc.cramfile),
        fna=lambda wc: f"zenodo/{'human' if 'human' in wc.cramfile else 'triticum'}.fna",
    output:
        bam="bamfiles/{cramfile}.bam",
        bai="bamfiles/{cramfile}.bam.bai",
    threads: 10
    shell:
        "samtools view -@ {threads} -T {input.fna} -b -o {output.bam} {input.cram} && "
        "samtools index -@ {threads} {output.bam}"
