### Dependencies

1. Install `rustc` >= 1.73
1. Install `matplotlib`, `Snakemake` and other optional dependencies you may need or want. For example: `conda create -n snkmk matplotlib snakemake snakemake-executor-plugin-slurm`

#### Download data

1. **Download the data from [here](https://zenodo.org/record/14760356)** (`10.5281/zenodo.14760356`), put these BAM files under `zenodo/` directory.
1. **Get GTF files from Ensembl for human and/ or wheat**. Put these under `regions/` directory and name them like we did.

This is how it should look like:

```
.
├── regions
│   ├── homo.v91.full.gtf
│   └── triticum.v60.full.gtf
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

> [!NOTE]
>
> There's also some extra data files in our zenodo repo, these were part of the paper but not the benchmark.

## Run benchmark

`snakemake --use-conda --benchmark-extended --config organism=homo ntimes=3 nthreads=8`

If you wanted to run this on an HPC cluster, we do provide an executor config file, just `$ run.sh` ;)

Optionally, you may want to adjust the different bin sizes at the top of the Snakefile, or even the hardcoded number of transcripts to be taken from the GTF.

> [!WARNING]
>
> The executor config file (`snk-slurm-exe/config.yaml`) doesn't specify number of CPUs. Instead, we are relying on the `Nthreads` parameter from `Snakefile`, same as when running locally. But you can adjust Slurm partition, memory, and runtime from there.

### Read results

> [!NOTE]
>
> Everything parsed from the results is saved as `*.csv` files.

Aside from `*walltime.png` and `*memory.png` boxplots (included in `output/report.html` on each run), there're other plots that may provide further insight into performance gains. Furthermore, there is a script to aggregate results from different configurations (`Nthreads`), you may use it as this:

```{bash}
./aggregate_runs.py --folders ~/data/bioinfo/benchmarking/hpc8_homo:8,~/data/bioinfo/benchmarking/hpc12_homo:12,~/data/bioinfo/benchmarking/hpc16_homo:16 --output memory_comparison --verbose
pandoc memory_comparison/memory_comparison_summary.md -o memory_comparison/summary.html
```
