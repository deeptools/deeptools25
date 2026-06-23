# deeptools 25 repo

## Setup
You need to have [rust](https://rustup.rs/) installed on your system before continuing, additionally, you need [pixi](https://pixi.sh/latest/) up and running too.

## Benchmark

Code for the benchmark is organized under the benchmark directory.
How to run:

## Benchmark

Run the benchmark with pixi:

  > To be implemented

## Example

To reproduce the example figures, [public data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183556) is sourced from [Sun et al.](https://doi.org/10.1038/s41586-023-06781-3).
Included are different ChIP data (H3K27ac, H3K4me1, H3K4me3, H3K9me3, H3k27me3), ATAC data, bisulfite sequencing data (DNA methylation) and RNA-seq data, in wild type NPC cells, or MSL2-KO NPC cells.

Recreate the example:

  > pixi run snakemake -s example/example.smk --configfile example/conf/smk_config.yaml --cores 10 -d /path/to/working/directory
