import pandas as pd
import numpy as np

# Down genes
downg = [] 
with open(snakemake.input.down, 'r') as f:
    for line in f:
        downg.append(line.strip().split('\t')[0])
downg = set(downg)

# Up genes
upg = [] 
with open(snakemake.input.up, 'r') as f:
    for line in f:
        upg.append(line.strip().split('\t')[0])
upg = set(upg)

# nonde genes
nondeg = []
with open(snakemake.input.nonde, 'r') as f:
    for line in f:
        nondeg.append(line.strip().split('\t')[0])
nondeg = set(nondeg)

# output files
downo = open(snakemake.output.downgtf, 'w')
downtss = open(snakemake.output.downbed, 'w')
upo = open(snakemake.output.upgtf, 'w')
uptss = open(snakemake.output.upbed, 'w')
nono = open(snakemake.output.nongtf, 'w')
nontss = open(snakemake.output.nonbed, 'w')

with open(snakemake.input.gtf, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        gid = None
        _l = line.strip().split('\t')[8].split(';')
        _li = line.strip().split('\t')
        for i in _l:
            if 'gene_id' in i:
                gid = i.split(' ')[1].replace('"', '')
        if gid in downg:
            downo.write(f'{line.strip()}\n')
            tss = None
            if _li[2] == 'gene':
                if _li[6] == '+':
                    tss = int(_li[3])
                    chrom = _li[0]
                else:
                    assert _li[6] == '-'
                    tss = int(_li[4])
                    chrom = _li[0]
                assert tss
                downtss.write(f'{chrom}\t{tss}\t{tss+1}\n')
        if gid in upg:
            upo.write(f'{line.strip()}\n')
            tss = None
            if _li[2] == 'gene':
                if _li[6] == '+':
                    tss = int(_li[3])
                    chrom = _li[0]
                else:
                    assert _li[6] == '-'
                    tss = int(_li[4])
                    chrom = _li[0]
                assert tss
                uptss.write(f'{chrom}\t{tss}\t{tss+1}\n')
        if gid in nondeg:
            nono.write(f'{line.strip()}\n')
            tss = None
            if _li[2] == 'gene':
                if _li[6] == '+':
                    tss = int(_li[3])
                    chrom = _li[0]
                else:
                    assert _li[6] == '-'
                    tss = int(_li[4])
                    chrom = _li[0]
                assert tss
                nontss.write(f'{chrom}\t{tss}\t{tss+1}\n')

downo.close()
downtss.close()
upo.close()
uptss.close()
nono.close()
nontss.close()

# For non-TSS ChIP data. We first identify peaks associated to downregulated gene, 
# If there are multiple, we retain the one with the largest delta between WT and KO.
def parse_chip(downg, of, uropdf, mbs, chip):
    of = open(of, 'w')
    df = pd.read_table(uropdf, low_memory=False)
    mat = np.load(mbs)
    matlabels = [bytes(x).decode('utf-8').rstrip('\x00').split('_')[2] for x in mat['labels']]
    ctrl_mask = np.array(matlabels) == 'ctrl'
    ko_mask = np.array(matlabels) == 'ko'

    for i in downg:
        tdf = df[df['gene_id'] == i]
        if len(tdf) == 0:
            continue
        else:
            region_indices = tdf.index.values
            matrix = mat['matrix'][region_indices] + 1
            ctrlm = np.median(matrix[:, ctrl_mask], axis=1)
            kom = np.median(matrix[:, ko_mask], axis=1)
            log2fc = np.log2(ctrlm / kom)
            ai = np.argmin(log2fc)
            if log2fc[ai] < snakemake.params.l2fc and chip in ['H3K27me3', 'H3K9me3']:
                ix = region_indices[ai]
                of.write(
                    f"{df['peak_chr'].iloc[ix]}\t{df['peak_start'].iloc[ix]}\t{df['peak_end'].iloc[ix]}\n"
                )
            if log2fc[ai] > snakemake.params.l2fc and chip in ['H3K4me1']:
                ix = region_indices[ai]
                of.write(
                    f"{df['peak_chr'].iloc[ix]}\t{df['peak_start'].iloc[ix]}\t{df['peak_end'].iloc[ix]}\n"
                )
    of.close()

parse_chip(downg, snakemake.output.k27_down, snakemake.input.uro_k27, snakemake.input.mbs_k27, 'H3K27me3')
parse_chip(downg, snakemake.output.k9_down, snakemake.input.uro_k9, snakemake.input.mbs_k9, 'H3K9me3')
parse_chip(downg, snakemake.output.k4_down, snakemake.input.uro_k4me1, snakemake.input.mbs_k4me1, 'H3K4me1')