[![DOI:10.5281/zenodo.21529137](https://img.shields.io/badge/DOI-10.5281/zenodo.21529137-yellow.svg)](https://doi.org/10.5281/zenodo.21529137)

# deeptools 25 repo

## Setup

You need to have [rust](https://rustup.rs/) installed on your system before continuing, additionally, you need [pixi](https://pixi.sh/latest/) and [conda](https://docs.conda.io/en/latest/) up and running too. Finally, your system needs internet access to download the required data.

### Benchmark

Run the benchmark:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda

Running specific functions only in the benchmark can be done by setting `what` in the config:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda --config what=bamCoverage

What can take values 'alignmentSieve', 'bamCoverage', 'bamCompare', 'multibamSummary' or 'computeMatrix'. If any other value is provided (or none is set), all benchmarks will be ran. The os flag will be set automatically, but can also be provided manually (used in final benchmark table to keep track of platform benchmark was run on):

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda --config what=bamCoverage os=my_fast_computer

An exhaustive run can be performed by setting `exhaustive` to true in the config:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda --config exhaustive=true


### Example

  > pixi run snakemake -s example/example.smk --cores 20 -d /path/to/working/directory --use-conda

## SLURM submission

In case you want to run either the benchmark of the example via slurm, an example profile is included in this repository. Appending `--profile slurm-example-config` results in jobs being submitted via SLURM rather than being executed locally. Make sure that the settings (`queue`, `mem_mb`, `runtime`, ...) in `slurm-example-config/config.yaml` are appropriate for your cluster.

## Collate figures

The results generated from both workflows (benchmark and example workflow), as well as the code to generate figure 1 are included in this repository. To collate/generate the figures in the `figures` directory, simply run:

  > pixi run bash figures/collate_figures.sh

## Public data
Both the benchmark as well as the example use publicly available data. The processed version of this data is included in the zenodo archive (and will be downloaded automatically by the respective workflows), the original source data can be retrieved through the accessions listed here.

| Purpose | Modality | Organism | SRA | BioProject | GEO | doi |
| --- | --- | --- | --- | --- | --- | --- |
| benchmark | ChIP - H3K9me3 | Human | SRR28592124 | PRJNA1097546 | GSE263436 | 10.1073/pnas.2412258121 |
| benchmark | RNA | Human | SRR28012902 | PRJNA1078228 | GSE256091 | 10.1158/0008-5472.CAN-23-1551 |
| benchmark | WGS | Human | SRR15494527 | PRJNA752223 | / | / |
| benchmark | ChIP - CENH3 | Wheat | SRR1686799 | PRJNA268976 | GSE63752 | 10.1105/tpc.19.00133 |
| benchmark | RNA | Wheat | SRR27822150 | PRJNA1071712 | GSE254797 | 10.3389/fpls.2024.1401135 |
| benchmark | WGS | Wheat | SRR27887047 | PRJNA956839 | / | 10.1038/s41586-024-07808-z |
| example | ChIP - H3K27ac wt1 | Mouse | SRR15829033 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K27ac wt2| Mouse | SRR15829034 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K27me3 wt1| Mouse | SRR15829036 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K27me3 wt2| Mouse | SRR15829037 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me1 wt1| Mouse | SRR15829042 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me1 wt2| Mouse | SRR15829043 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me3 wt1| Mouse | SRR15829045 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me3 wt2| Mouse | SRR15829046 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K9me3 wt1| Mouse | SRR15829048 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K9me3 wt2| Mouse | SRR15829049 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - input wt1| Mouse | SRR15829058 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - input wt2| Mouse | SRR15829059 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | RNA wt1| Mouse | SRR15829077 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | RNA wt2| Mouse | SRR15829078 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ATAC wt1| Mouse | SRR15829121 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ATAC wt2| Mouse | SRR15829122 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | WGBS wt1| Mouse | SRR15829019 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | WGBS wt2| Mouse | SRR15829020 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | RNA ko1| Mouse | SRR15829149 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | RNA ko2| Mouse | SRR15829150 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ATAC wt1| Mouse | SRR15829084 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ATAC wt2| Mouse | SRR15829085 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | WGBS wt1| Mouse | SRR15829089 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | WGBS wt2| Mouse | SRR15829090 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K27ac ko1| Mouse | SRR15829101 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K27ac ko2| Mouse | SRR15829102 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K27me3 ko1| Mouse | SRR15829103 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K27me3 ko2| Mouse | SRR15829104 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me1 ko1| Mouse | SRR15829107 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me1 ko2| Mouse | SRR15829108 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me3 ko1| Mouse | SRR15829109  | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K4me3 ko2| Mouse | SRR15829110 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K9me3 ko1| Mouse | SRR15829111 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - H3K9me3 ko2| Mouse | SRR15829124 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - input ko1| Mouse | SRR15829131 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
| example | ChIP - input ko2| Mouse | SRR15829132 | PRJNA762170 | GSE183764 | 10.1038/s41586-023-06781-3 |
