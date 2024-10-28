# deeptools25

> repo for benchmarking, comparison, test data and reference construction for the 25 remake.

## Install

```{bash}
sudo apt install clang
conda create -n deeptools25 python==3.9.12 deeptools==3.5.5 pybigwig==0.3.22 snakemake==8.24.1
git clone git@github.com:WardDeb/deepTools.git ../deepToolsWard && cd ../deepToolsWard
git checkout -b maturin && pip install -e .
```

### Run benchmark

```{bash}
cp -r /data/manke/processing/deboutte/tmp/testfiles /tmp/  # TODO: find public test files (e.g. from bigtools paper)
snakemake --forcerun --snakefile scripts/bench.snk --benchmark-extended --cores 1
rm ./output/*.{png,txt} && cp -v /tmp/testfiles/output/*.{png,txt} ./output/
```

> [!NOTE]
> Adjust `./bench.cfg` as desired.
