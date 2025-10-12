# deeptools 25 repo

## Setup
You need to have [rust](https://rustup.rs/) installed on your system before continuing, as well as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

From there, you can set up the conda environment

 > conda env create -f conf/env.yml -n deeptools_benchmark  
 > conda activate deeptools_benchmark  
 > CFLAGS=-I/$CONDA_PREFIX/include LIBCLANG_PATH=$CONDA_PREFIX/lib pip install -r conf/requirements.txt

Alternatively, if you have [pixi](https://pixi.sh/latest/) available, a pixi configuration file is provided too.

## Benchmark

Code for the benchmark is organized under the benchmark directory.
How to run:

## Example

To reproduce the example figures, [public data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183556) is sourced from [Sun et al.](https://doi.org/10.1038/s41586-023-06781-3).
Included are different ChIP data (H3K27ac, H3K4me1, H3K4me3, H3K9me3, H3k27me3), ATAC data, bisulfite sequencing data (DNA methylation) and RNA-seq data, in wild type NPC cells, or MSL2-KO NPC cells.

### Set up

There are two entry points to reproduce the figures:
  - from the raw data  
  - from CRAM files in zenodo  

Assumed is you have the `deeptools_benchmark` conda environment installed (see above).
In case you want to start from the raw data, the workflow here assumes you have a working [snakePipes (>= 3.1.0)](https://github.com/maxplanck-ie/snakepipes) version installed and working, in a specific conda environment. The settings can be set in the `example/conf/smk_config.yaml` file. Note that both modes require compute node access to the internet. By default the mode is set to `zenodo`, and requires no additional parameters to be set.

### Running
With the `deeptools_benchmark` conda environment activated, you can reproduce the figures:

with conda environment (deeptools_benchmark) activated:

  > snakemake -s example/example.smk --configfile example/conf/smk_config.yaml --cores 10 -d /path/to/working/directory

Or with pixi:

  > pixi run snakemake -s example/example.smk --configfile example/conf/smk_config.yaml --cores 10 -d /path/to/working/directory

