import numpy as np
import pandas as pd

# Down genes
downg = []
with open(snakemake.input.down, "r") as f:
    for line in f:
        downg.append(line.strip().split("\t")[0])
downg = set(downg)

# Up genes
upg = []
with open(snakemake.input.up, "r") as f:
    for line in f:
        upg.append(line.strip().split("\t")[0])
upg = set(upg)

# nonde genes
nondeg = []
with open(snakemake.input.nonde, "r") as f:
    for line in f:
        nondeg.append(line.strip().split("\t")[0])
nondeg = set(nondeg)


def parse_chip(deg, ofs, uropdf, mbs, chip):
    downg, upg, nondeg = deg
    downo, upo, nondeo = ofs
    downo = open(downo, "w")
    upo = open(upo, "w")
    nondeo = open(nondeo, "w")

    df = pd.read_table(uropdf, low_memory=False)
    mat = np.load(mbs)
    matlabels = [
        bytes(x).decode("utf-8").rstrip("\x00").split("_")[2] for x in mat["labels"]
    ]
    ctrl_mask = np.array(matlabels) == "ctrl"
    ko_mask = np.array(matlabels) == "ko"

    for i in downg:
        tdf = df[df["gene_id"] == i]
        if len(tdf) == 0:
            continue
        else:
            region_indices = tdf.index.values
            matrix = mat["matrix"][region_indices] + 1
            ctrlm = np.median(matrix[:, ctrl_mask], axis=1)
            kom = np.median(matrix[:, ko_mask], axis=1)
            log2fc = np.log2(ctrlm / kom)
            if chip in ["H3K27me3", "H3K9me3"]:
                ai = np.argmin(log2fc)
                if log2fc[ai] <= -snakemake.params.l2fc:
                    ix = region_indices[ai]
                    downo.write(
                        f"{df['peak_chr'].iloc[ix]}\t{df['peak_start'].iloc[ix]}\t{df['peak_end'].iloc[ix]}\n"
                    )
            else:
                ai = np.argmax(log2fc)
                if log2fc[ai] >= snakemake.params.l2fc:
                    ix = region_indices[ai]
                    downo.write(
                        f"{df['peak_chr'].iloc[ix]}\t{df['peak_start'].iloc[ix]}\t{df['peak_end'].iloc[ix]}\n"
                    )

    for i in upg:
        tdf = df[df["gene_id"] == i]
        if len(tdf) == 0:
            continue
        else:
            region_indices = tdf.index.values
            matrix = mat["matrix"][region_indices] + 1
            ctrlm = np.median(matrix[:, ctrl_mask], axis=1)
            kom = np.median(matrix[:, ko_mask], axis=1)
            log2fc = np.log2(ctrlm / kom)
            if chip in ["H3K27me3", "H3K9me3"]:
                ai = np.argmax(log2fc)
                if log2fc[ai] >= snakemake.params.l2fc:
                    ix = region_indices[ai]
                    upo.write(
                        f"{df['peak_chr'].iloc[ix]}\t{df['peak_start'].iloc[ix]}\t{df['peak_end'].iloc[ix]}\n"
                    )
            else:
                ai = np.argmin(log2fc)
                if log2fc[ai] <= -snakemake.params.l2fc:
                    ix = region_indices[ai]
                    upo.write(
                        f"{df['peak_chr'].iloc[ix]}\t{df['peak_start'].iloc[ix]}\t{df['peak_end'].iloc[ix]}\n"
                    )

    # Genes nonde
    for i in nondeg:
        tdf = df[df["gene_id"] == i]
        if len(tdf) == 0:
            continue
        else:
            region_indices = tdf.index.values
            matrix = mat["matrix"][region_indices] + 1
            ctrlm = np.median(matrix[:, ctrl_mask], axis=1)
            kom = np.median(matrix[:, ko_mask], axis=1)
            log2fc = np.log2(kom / ctrlm)

            ai = np.argmin(np.abs(log2fc))
            if np.abs(log2fc[ai]) < snakemake.params.l2fc:
                ix = region_indices[ai]
                nondeo.write(
                    f"{df['peak_chr'].iloc[ix]}\t{df['peak_start'].iloc[ix]}\t{df['peak_end'].iloc[ix]}\n"
                )
    downo.close()
    upo.close()
    nondeo.close()


parse_chip(
    (downg, upg, nondeg),
    (snakemake.output.down, snakemake.output.up, snakemake.output.nonde),
    snakemake.input.urop,
    snakemake.input.mbs,
    snakemake.params.chip,
)
