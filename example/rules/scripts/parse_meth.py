import pandas as pd

downg = set(l.strip().split('\t')[0] for l in open(snakemake.input.down))
upg = set(l.strip().split('\t')[0] for l in open(snakemake.input.up))
nondeg = set(l.strip().split('\t')[0] for l in open(snakemake.input.nonde))

bins = pd.read_table(snakemake.input.bins, header=None, names=['chrom', 'start', 'end', 'gene_id'])
counts = pd.read_table(snakemake.input.counts)
counts.columns = [c.strip("'#") for c in counts.columns]
assert len(bins) == len(counts)

samplecols = [c for c in counts.columns if c not in ('chr', 'start', 'end')]
ctrl_cols = [c for c in samplecols if '_ctrl_' in c]
ko_cols = [c for c in samplecols if '_ko_' in c]

bins['diff'] = counts[ko_cols].mean(axis=1).values - counts[ctrl_cols].mean(axis=1).values
mdiff = snakemake.params.mdiff


def strongest_bin(genes, want_extreme):
    rows = []
    for g in genes:
        sub = bins[bins['gene_id'] == g]
        if sub.empty:
            continue
        ix = sub['diff'].abs().idxmax() if want_extreme else sub['diff'].abs().idxmin()
        passes = abs(sub.loc[ix, 'diff']) >= mdiff if want_extreme else abs(sub.loc[ix, 'diff']) < mdiff
        if passes:
            rows.append(sub.loc[ix])
    return rows


def write_bed(path, rows):
    with open(path, 'w') as f:
        for row in rows:
            f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\n")


write_bed(snakemake.output.down, strongest_bin(downg, True))
write_bed(snakemake.output.up, strongest_bin(upg, True))
write_bed(snakemake.output.nonde, strongest_bin(nondeg, False))
