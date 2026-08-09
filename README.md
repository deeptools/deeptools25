# deeptools 25 repo

## Setup

You need to have [rust](https://rustup.rs/) installed on your system before continuing, additionally, you need [pixi](https://pixi.sh/latest/) and [conda](https://docs.conda.io/en/latest/) up and running too. Finally, your system needs internet access to download the required data.

Run the benchmark:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda

Running specific functions only in the benchmark can be done by setting `what` in the config:

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda --config what=bamCoverage

What can take values 'alignmentSieve', 'bamCoverage', 'bamCompare', 'multibamSummary' or 'computeMatrix'. If any other value is provided (or none is set), all benchmarks will be ran. The os flag will be set automatically, but can also be provided manually (used in final benchmark table to keep track of platform benchmark was run on):

  > pixi run snakemake -s benchmark/benchmark.smk --cores 20 -d /path/to/working/directory --use-conda --config what=bamCoverage os=my_fast_computer

Run the example:

  > pixi run snakemake -s example/example.smk --cores 20 -d /path/to/working/directory --use-conda

## Running on SLURM

First, adjust `cpus_per_task` to the number of cores of the biggest node on the cluster (or partition, but that would require extra config.) Then, similar to above, plus `--profile benchmark/conf/slurm-profile`:

  > pixi run snakemake -s benchmark/benchmark.smk -d /path/to/working/directory --profile benchmark/conf/slurm-profile

The profile adds the SLURM executor, and reserves a whole node per job (combining `--exclusive` + `cpus_per_task`), and doesn't cap memory further since the node isn't shared with anyone else. (The rules' own `resources:` -e.g. `mem_mb`- stay as-is and still apply.)

Repeats default to 3 (set in `benchmark/benchmark.smk`); the profile bumps that to 10 via `config: [repeats=10]`, since SLURM has the room and queue time to afford more replicates per tool/version/sample. Override on the command line with `--config repeats=N` if you want something else.
