# deepTools benchmark

## Setup
You need to have [rust](https://rustup.rs/) installed on your system before continuing, as well as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

From there, you can set up the conda environment

 > conda env create -f conf/env.yml -n deeptools_benchmark  
 > conda activate deeptools_benchmark
 > CFLAGS=-I/$CONDA_PREFIX/include LIBCLANG_PATH=$CONDA_PREFIX/lib pip install -r conf/requirements.txt

Subsequently, run the benchmark.

 > snakemake -d working_directory -s path_to_repo/benchmark.smk --use-conda --config os=RHEL

You could either specify the number of cores (for a local run, with --cores 40) or a profile for distributed run. Note that --use-conda is required. The --config os is optional, but is included as a column in the benchmark output csv. This is useful to compare different OS'es. If os=osx, then the snakemake LOG file will be parsed to get the used memory (as this doesn't work with snakemake benchmark atm).
