# deeptools 25 repo

## Setup
You need to have [rust](https://rustup.rs/) installed on your system before continuing, additionally, you need [pixi](https://pixi.sh/latest/) and [conda](https://docs.conda.io/en/latest/) up and running too. Finally, your system needs internet access to download the required data.

Run the benchmark:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda

Running specific functions only in the benchmark can be done by setting `what` in the config:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda --config what=bamCoverage

What can take values 'bamCoverage', 'bamCompare', 'multibamSummary' or 'computeMatrix'. If any other value is provided (or none is set), all benchmarks will be ran. The os flag will be set automatically, but can also be provided manually (used in final benchmark table to keep track of platform benchmark was run on):

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda --config what=bamCoverage os=my_fast_computer

Run the example:

  > pixi run snakemake -s example/example.smk --cores 20 -d /path/to/working/directory --use-conda
