### Preparation

```{bash}
conda create -n snakemake snakemake
conda activate snakemake
```

#### Download data

> WIP.

## Run benchmark

> [!NOTE]
> First, adjust `./bench.cfg` as desired with a symlink to any of the CFG files bundled within,
> for example, to use our `olddata/`: `ln -sf old.cfg bench.cfg`

`snakemake --forcerun --use-conda --benchmark-extended --cores 12`


#### Clean

```{bash}
rm ./output/*.{png,txt}
cp -v /tmp/testfiles/output/*.{png,txt} ./output/
```
