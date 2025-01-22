# deeptools25

> repo for benchmarking, comparison, test data and reference construction for the 25 remake.

## Install

```{bash}
mamba create -n deeptools25 clang-19 python~=3.9 deeptools~=3.5 pybigwig~=0.3 snakemake~=8.24
git clone git@github.com:WardDeb/deepTools.git deepToolsWard
cd deepToolsWard && git checkout -b maturin && pip install -e .
```

### Run benchmark

```{bash}
cp -r /data/manke/processing/deboutte/tmp/testfiles /tmp/
snakemake --forcerun --use-conda --snakefile scripts/bench.snk --benchmark-extended --cores 1
rm ./output/*.{png,txt} && cp -v /tmp/testfiles/output/*.{png,txt} ./output/
```

> [!NOTE]
> Adjust `./bench.cfg` as desired (E.g. a symlink to any of the CFG files bundled within)
