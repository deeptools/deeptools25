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

The profile adds the SLURM executor and reserves a whole node per job (combining `--exclusive` + `cpus_per_task`). It also uses `set-threads:` to make every tool actually use that many threads (`-p {threads}`/`-@ {threads}`) -- without it, each rule would keep using whatever thread count is hardcoded in the Snakefile (often much lower), leaving most of the reserved node idle.

Two profile variants exist, differing only in how many threads each tool gets while still reserving the same exclusive node:
- `benchmark/conf/slurm-profile` -- every rule uses the full node core count (e.g. 72 threads on a 72-core node). Maximizes throughput per job, but isn't representative of how most users actually submit jobs.
- `benchmark/conf/slurm-profile-24t` -- caps each tool at 24 threads, closer to a typical user's request, while still landing on the same exclusive/big-memory node. Use this for a more realistic HPC measurement.

A whole node not being shared with anyone else means memory doesn't need capping to protect other tenants -- but SLURM still enforces `mem_mb` as a hard per-job cgroup limit regardless of node size, and cranking up thread count inflates a tool's per-thread buffer/chunk memory well past what the Snakefile's default `resources:` assume (those were sized for low thread counts). Both profiles' `set-resources:` blocks bump `mem_mb` for the specific rules we've hit this on in practice (`bamcompare_dt4`, `multibamsummary_dt3`, `computeMatrix_dt4` on the largest input); if you push thread counts further, watch for `oom_kill` in the SLURM job logs and raise `mem_mb` for whichever rule hits it.

Repeats default to 3 (set in `benchmark/benchmark.smk`); the profiles bump that to 10 via `config: [repeats=10]`, since SLURM has the room and queue time to afford more replicates per tool/version/sample. Override on the command line with `--config repeats=N` if you want something else.
