import numpy as np
import pandas as pd

mode = snakemake.params.mode
tool = snakemake.params.tool
sample = snakemake.params.sample

row = {"tool": tool, "sample": sample}

if mode == "bigwig":
    import pyBigWig

    bw3 = pyBigWig.open(snakemake.input.dt3)
    bw4 = pyBigWig.open(snakemake.input.dt4)
    v3, v4 = [], []
    for chrom, length in bw3.chroms().items():
        if chrom not in bw4.chroms():
            continue
        v3.append(np.nan_to_num(np.array(bw3.values(chrom, 0, length))))
        v4.append(np.nan_to_num(np.array(bw4.values(chrom, 0, length))))
    bw3.close()
    bw4.close()
    v3 = np.concatenate(v3)
    v4 = np.concatenate(v4)
    row["n_values"] = len(v3)
    row["pearson_r"] = float(np.corrcoef(v3, v4)[0, 1])
    row["max_abs_diff"] = float(np.max(np.abs(v3 - v4)))

elif mode == "npz":
    d3 = np.load(snakemake.input.dt3)["matrix"]
    d4 = np.load(snakemake.input.dt4)["matrix"]
    row["n_values"] = d3.size
    row["pearson_r"] = float(np.corrcoef(d3.flatten(), d4.flatten())[0, 1])
    row["max_abs_diff"] = float(np.max(np.abs(d3 - d4)))

elif mode == "matrix":
    from deeptools.heatmapper import heatmapper

    hm3 = heatmapper()
    hm3.read_matrix_file(snakemake.input.dt3)
    hm4 = heatmapper()
    hm4.read_matrix_file(snakemake.input.dt4)
    m3 = np.ma.filled(hm3.matrix.matrix, 0.0)
    m4 = np.ma.filled(hm4.matrix.matrix, 0.0)
    row["n_values"] = m3.size
    row["pearson_r"] = float(np.corrcoef(m3.flatten(), m4.flatten())[0, 1])
    row["max_abs_diff"] = float(np.max(np.abs(m3 - m4)))

elif mode == "bam":
    import pysam

    def read_index(path):
        idx = {}
        with pysam.AlignmentFile(path, "rb") as f:
            for read in f:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                idx[(read.query_name, read.is_read1)] = (
                    read.reference_start,
                    read.query_length,
                    read.flag,
                )
        return idx

    idx3 = read_index(snakemake.input.dt3)
    idx4 = read_index(snakemake.input.dt4)
    common = idx3.keys() & idx4.keys()
    row["n_reads_dt3"] = len(idx3)
    row["n_reads_dt4"] = len(idx4)
    row["n_common_reads"] = len(common)
    row["n_pos_mismatch"] = sum(idx3[k][0] != idx4[k][0] for k in common)
    row["n_qlen_mismatch"] = sum(idx3[k][1] != idx4[k][1] for k in common)
    row["n_flag_mismatch"] = sum(idx3[k][2] != idx4[k][2] for k in common)

else:
    raise ValueError(f"unknown comparison mode: {mode}")

pd.DataFrame([row]).to_csv(snakemake.output.tsv, sep="\t", index=False)
