# deeptools 25 repo

## Setup
You need to have [rust](https://rustup.rs/) installed on your system before continuing, additionally, you need [pixi](https://pixi.sh/latest/) and [conda](https://docs.conda.io/en/latest/) up and running too. Finally, your system needs internet access to download the required data.

Run the benchmark:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda

Run the example:

  > pixi run snakemake -s example/example.smk --cores 20 -d /path/to/working/directory --use-conda
