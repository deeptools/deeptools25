### Dependencies

1. Install `rustc` >= 1.73
1. Install `Clang`, `matplotlib` and `Snakemake`. For example: `conda create -n snkmk clang-19 matplotlib snakemake bedtools ucsc-bigwiginfo`

> [!NOTE]
> Alternatively, you can use any modern IDE with [Devcontainers](https://code.visualstudio.com/docs/devcontainers/containers) support.

#### Download data

1. **Download the data from [here](https://zenodo.org/record/14760356)** (`10.5281/zenodo.14760356`), put these BAM files under `zenodo/` directory.
1. **Get GTF files from Ensembl for human and/ or wheat**. Put these under `regions/` directory. If you are not inclined into benchmarking on full transcriptome, downsample those files like this: `grep 'transcript_id' homo.v91.full.gtf | shuf | head -n 25000 | bedtools sort -i - > homo.v91.sample25k.gtf`

This is how it should look like:

```
.
├── regions
│   ├── homo.v91.full.gtf
│   ├── homo.v91.sample25k.gtf
│   ├── triticum.v60.full.gtf
│   └── triticum.v60.sample25k.gtf
└── zenodo
    ├── bigwigs
    │   ├── human_chip_SRR28592124.bw
    │   ├── human_chip_SRR28592125.bw
    │   ├── human_chip_SRR28592131.bw
    │   ├── human_chip_SRR28592132.bw
    │   ├── human_rna_SRR28012902.bw
    │   ├── human_rna_SRR28012903.bw
    │   ├── human_rna_SRR28012904.bw
    │   ├── human_rna_SRR28012905.bw
    │   └── human_wgs_SRR15494527.bw
    ├── human_chip_SRR28592124.bam
    ├── human_chip_SRR28592124.bam.bai
    ├── human_rna_SRR28012902.bam
    ├── human_rna_SRR28012902.bam.bai
    ├── human_wgs_SRR15494527.bam
    ├── human_wgs_SRR15494527.bam.bai
    ├── triticum_chip_SRR1686799.bam
    ├── triticum_chip_SRR1686799.bam.csi
    ├── triticum_rna_SRR27822150.bam
    ├── triticum_rna_SRR27822150.bam.csi
    ├── triticum_wgs_SRR27887047.bam
    └── triticum_wgs_SRR27887047.bam.csi
```

## Run benchmark

Adjust the variables at the top of `Snakefile` however you like it, and then just `$ run.sh` ;)
