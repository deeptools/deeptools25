# deeptools 25 repo

## Setup
You need to have [rust](https://rustup.rs/) installed on your system before continuing, as well as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

From there, you can set up the conda environment

 > conda env create -f conf/env.yml -n deeptools_benchmark  
 > conda activate deeptools_benchmark  
 > CFLAGS=-I/$CONDA_PREFIX/include LIBCLANG_PATH=$CONDA_PREFIX/lib pip install -r conf/requirements.txt

## Benchmark

Code for the benchmark is organized under the benchmark directory.

## Example

Example run for deepTools usage on multimodal data is available under the example directory.
