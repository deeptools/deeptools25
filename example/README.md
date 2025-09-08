# Example usage

To reproduce the example figure, [public data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183556) is sourced from [Sun et al.](https://doi.org/10.1038/s41586-023-06781-3).
Included are different ChIP data (H3K27ac, H3K4me1, H3K4me3, H3K9me3, H3k27me3), ATAC data, bisulfite sequencing data (DNA methylation) and RNA-seq data, in wild type NPC cells, or MSL2-KO NPC cells.

## Set up

There are two entry points to reproduce the figures:
  - from the raw data  
  - from CRAM files in zenodo  

Assumed is you have the `deeptools25` conda environment installed. 
In case you want to start from the raw data, the workflow here assumes you have a working [snakePipes (>= 3.1.0)](https://github.com/maxplanck-ie/snakepipes) version installed and working, in a specific conda environment. The settings can be set in the `conf/smk_config.yaml` file. Note that both modes require compute node access to the internet.

## Running

snakemake -s example.smk --configfile conf/smk_config.yaml --cores 10 -d /path/to/working/directory --use-conda 