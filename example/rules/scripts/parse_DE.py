import pandas as pd
upg = []
downg = [] 
with open(snakemake.input.up, 'r') as f:
    for line in f:
        upg.append(line.strip().split('\t')[0])
with open(snakemake.input.down, 'r') as f:
    for line in f:
        downg.append(line.strip().split('\t')[0])

upg = set(upg)
downg = set(downg)

upo = open(snakemake.output.upgtf, 'w')
downo = open(snakemake.output.downgtf, 'w')
uptss = open(snakemake.output.upbed, 'w')
downtss = open(snakemake.output.downbed, 'w')

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

upo.close()
downo.close()
uptss.close()
downtss.close()

k27uo = open(snakemake.output.k27_up, 'w')
k9uo = open(snakemake.output.k9_up, 'w')
k27do = open(snakemake.output.k27_down, 'w')
k9do = open(snakemake.output.k9_down, 'w')

k27df = pd.read_table(snakemake.input.uro_k27, low_memory=False)
for i in upg:
    tdf = k27df[k27df['gene_id'] == i]
    if len(tdf) == 0:
        continue
    if len(tdf) == 1:
        k27uo.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )
    else:
        tdf = tdf.loc[[tdf['distance'].idxmin()]]
        k27uo.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )

for i in downg:
    tdf = k27df[k27df['gene_id'] == i]
    if len(tdf) == 0:
        continue
    if len(tdf) == 1:
        k27do.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )
    else:
        tdf = tdf.loc[[tdf['distance'].idxmin()]]
        k27do.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )




k9df = pd.read_table(snakemake.input.uro_k9, low_memory=False)
for i in upg:
    tdf = k9df[k9df['gene_id'] == i]
    if len(tdf) == 0:
        continue
    if len(tdf) == 1:
        k9uo.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )
    else:
        tdf = tdf.loc[[tdf['distance'].idxmin()]]
        k9uo.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )

for i in downg:
    tdf = k9df[k9df['gene_id'] == i]
    if len(tdf) == 0:
        continue
    if len(tdf) == 1:
        k9do.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )
    else:
        tdf = tdf.loc[[tdf['distance'].idxmin()]]
        k9do.write(
            f"{tdf['peak_chr'].iloc[0]}\t{tdf['peak_start'].iloc[0]}\t{tdf['peak_end'].iloc[0]}\n"
        )
