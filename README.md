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

Example run for deepTools usage on multimodal data is available under the example directory.
Note that this workflow requires downloading data from either zenodo or SRA (default is zenodo).
If you prefer to start from SRA / raw fastq files, a working installation of snakePipes need to be available in a conda environment. The appropriate settings need to be filled out in `example/conf/smk_config.yaml`.

To regenerate the figures:

via conda:

  > conda activate deeptools_benchmark
  > snakemake -s example/example.smk --configfile example/conf/smk_config.yaml --cores 40 -d WORKDIR --use-conda

or with pixi:

  > pixi run snakemake -s example/example.smk --configfile example/conf/smk_config.yaml --cores 40 -d WORKDIR
