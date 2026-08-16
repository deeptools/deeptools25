import io
import subprocess
import pandas as pd

bam = snakemake.input.bam
sample = snakemake.wildcards.cramfile

def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout

flagstat = run(f"samtools flagstat -@ {snakemake.threads} {bam}")
total = int(flagstat.splitlines()[0].split()[0])
dup = int([l for l in flagstat.splitlines() if "duplicates" in l][0].split()[0])
mapped_pct = float(flagstat.split("mapped (")[1].split("%")[0])
paired_pct = float(flagstat.split("properly paired (")[1].split("%")[0])

stats = run(f"samtools stats -@ {snakemake.threads} {bam}")
sn = {l.split("\t")[1].rstrip(":"): l.split("\t")[2] for l in stats.splitlines() if l.startswith("SN")}

coverage = run(f"samtools coverage {bam}")
cov = pd.read_csv(io.StringIO(coverage), sep="\t")
contig_len = cov["endpos"] - cov["startpos"] + 1
genome_len = contig_len.sum()
genome_cov_pct = 100 * cov["covbases"].sum() / genome_len
mean_depth = (cov["meandepth"] * contig_len).sum() / genome_len

row = pd.DataFrame([{
    "sample": sample,
    "total_reads": total,
    "pct_mapped": mapped_pct,
    "pct_properly_paired": paired_pct,
    "pct_duplicates": 100 * dup / total if total else 0.0,
    "avg_read_length": float(sn["average length"]),
    "insert_size_mean": float(sn.get("insert size average", "nan")),
    "insert_size_sd": float(sn.get("insert size standard deviation", "nan")),
    "error_rate": float(sn["error rate"]),
    "genome_coverage_pct": genome_cov_pct,
    "mean_depth": mean_depth,
}])
row.to_csv(snakemake.output.tsv, sep="\t", index=False)
