import numpy as np
import pandas as pd

mode = snakemake.params.mode
tool = snakemake.params.tool
sample = snakemake.params.sample

row = {"tool": tool, "sample": sample}

CHUNK_SIZE = 50_000_000


def streaming_compare(chunk_pairs):
    """
    chunk_pairs: iterable of (x, y) same-length float32 numpy arrays.
    Returns (n, pearson_r, max_abs_diff) without ever materializing full arrays.
    """
    n = 0
    sx = sy = sxx = syy = sxy = 0.0
    max_abs_diff = 0.0
    for x, y in chunk_pairs:
        n += x.size
        sx += x.sum(dtype=np.float64)
        sy += y.sum(dtype=np.float64)
        sxx += (x * x).sum(dtype=np.float64)
        syy += (y * y).sum(dtype=np.float64)
        sxy += (x * y).sum(dtype=np.float64)
        chunk_max = float(np.abs(x - y).max())
        if chunk_max > max_abs_diff:
            max_abs_diff = chunk_max

    if n == 0:
        return 0, float("nan"), float("nan")

    cov = sxy / n - (sx / n) * (sy / n)
    var_x = sxx / n - (sx / n) ** 2
    var_y = syy / n - (sy / n) ** 2
    denom = np.sqrt(var_x * var_y)
    pearson_r = float(cov / denom) if denom > 0 else float("nan")
    return n, pearson_r, max_abs_diff


if mode == "bigwig":
    import pyBigWig

    def chrom_chunks(bw, chrom, length):
        for start in range(0, length, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, length)
            try:
                vals = bw.values(chrom, start, end, numpy=True)
            except TypeError:
                vals = np.array(bw.values(chrom, start, end))
            yield np.nan_to_num(vals, copy=False).astype(np.float32, copy=False)

    def paired_chunks(bw3, bw4):
        for chrom, length in bw3.chroms().items():
            if chrom not in bw4.chroms():
                continue
            yield from zip(chrom_chunks(bw3, chrom, length), chrom_chunks(bw4, chrom, length))

    bw3 = pyBigWig.open(snakemake.input.dt3)
    bw4 = pyBigWig.open(snakemake.input.dt4)
    n_values, pearson_r, max_abs_diff = streaming_compare(paired_chunks(bw3, bw4))
    bw3.close()
    bw4.close()

    row["n_values"] = n_values
    row["pearson_r"] = pearson_r
    row["max_abs_diff"] = max_abs_diff


elif mode == "npz":
    d3 = pd.read_csv(snakemake.input.dt3, sep="\t", header=0, low_memory=False)
    d4 = pd.read_csv(snakemake.input.dt4, sep="\t", header=0, low_memory=False)

    row["shape_dt3"] = str(d3.shape)
    row["shape_dt4"] = str(d4.shape)
    row["shape_mismatch"] = d3.shape != d4.shape

    def coord_keys(df):
        return list(zip(df.iloc[:, 0].astype(str),
                        df.iloc[:, 1].astype(np.int64),
                        df.iloc[:, 2].astype(np.int64)))

    pos4 = {c: j for j, c in enumerate(coord_keys(d4))}
    matches = [(i, pos4[c]) for i, c in enumerate(coord_keys(d3)) if c in pos4]

    idx = np.linspace(0, len(matches) - 1, min(len(matches), 20_000)).astype(int)
    sel = [matches[k] for k in idx]
    v3_rows = [i for i, _ in sel]
    v4_rows = [j for _, j in sel]

    v3 = d3.iloc[v3_rows, 3:].to_numpy(dtype=np.float64)
    v4 = d4.iloc[v4_rows, 3:].to_numpy(dtype=np.float64)

    row["n_bins_dt3"] = len(d3)
    row["n_bins_dt4"] = len(d4)
    row["n_common_bins"] = len(matches)
    row["n_common_bins_sampled"] = len(sel)
    row["n_values"] = v3.size
    row["pearson_r"] = float(np.corrcoef(v3.ravel(), v4.ravel())[0, 1]) if v3.size else float("nan")
    row["max_abs_diff"] = float(np.max(np.abs(v3 - v4))) if v3.size else float("nan")


elif mode == "matrix":
    from deeptools.heatmapper import heatmapper

    ROW_LIMIT = 20_000

    hm3 = heatmapper()
    hm3.read_matrix_file(snakemake.input.dt3)
    hm4 = heatmapper()
    hm4.read_matrix_file(snakemake.input.dt4)
    m3 = hm3.matrix.matrix
    m4 = hm4.matrix.matrix
    row["shape_dt3"] = str(m3.shape)
    row["shape_dt4"] = str(m4.shape)
    row["shape_mismatch"] = m3.shape != m4.shape

    n_rows = min(m3.shape[0], m4.shape[0], ROW_LIMIT)
    m3_head = np.ma.filled(m3[:n_rows], 0.0)
    m4_head = np.ma.filled(m4[:n_rows], 0.0)

    row["n_values"] = m3_head.size
    row["pearson_r"] = float(np.corrcoef(m3_head.ravel(), m4_head.ravel())[0, 1])
    row["max_abs_diff"] = float(np.max(np.abs(m3_head - m4_head)))

elif mode == "bam":
    import pysam

    INDEX_LIMIT = 20_000

    def sample_index(path, limit):
        idx = {}
        n_reads = 0
        with pysam.AlignmentFile(path, "rb") as f:
            for read in f:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                n_reads += 1
                if len(idx) < limit:
                    idx[(read.query_name, read.is_read1)] = (
                        read.reference_start,
                        read.query_length,
                        read.flag,
                    )
        return idx, n_reads

    def targeted_index(path, targets):
        idx = {}
        n_reads = 0
        with pysam.AlignmentFile(path, "rb") as f:
            for read in f:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                n_reads += 1
                key = (read.query_name, read.is_read1)
                if key in targets:
                    idx[key] = (read.reference_start, read.query_length, read.flag)
        return idx, n_reads

    idx3, n_reads_dt3 = sample_index(snakemake.input.dt3, INDEX_LIMIT)
    idx4, n_reads_dt4 = targeted_index(snakemake.input.dt4, idx3)

    common = idx3.keys() & idx4.keys()
    row["n_reads_dt3"] = n_reads_dt3
    row["n_reads_dt4"] = n_reads_dt4
    row["n_common_reads_sampled"] = len(common)
    row["n_pos_mismatch"] = sum(idx3[k][0] != idx4[k][0] for k in common)
    row["n_qlen_mismatch"] = sum(idx3[k][1] != idx4[k][1] for k in common)
    row["n_flag_mismatch"] = sum(idx3[k][2] != idx4[k][2] for k in common)


else:
    raise ValueError(f"unknown comparison mode: {mode}")

pd.DataFrame([row]).to_csv(snakemake.output.tsv, sep="\t", index=False)
