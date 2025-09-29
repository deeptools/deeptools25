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