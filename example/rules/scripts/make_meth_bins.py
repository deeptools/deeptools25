downg = set(l.strip().split('\t')[0] for l in open(snakemake.input.down))
upg = set(l.strip().split('\t')[0] for l in open(snakemake.input.up))
nondeg = set(l.strip().split('\t')[0] for l in open(snakemake.input.nonde))
allgenes = downg | upg | nondeg

binsize = snakemake.params.binsize
flank = snakemake.params.flank

with open(snakemake.input.gtf) as f, open(snakemake.output.bed, 'w') as out:
    for line in f:
        if line.startswith('#'):
            continue
        li = line.strip().split('\t')
        if li[2] != 'gene':
            continue
        gid = None
        for a in li[8].split(';'):
            if 'gene_id' in a:
                gid = a.strip().split(' ')[1].replace('"', '')
        if gid not in allgenes:
            continue
        chrom = li[0]
        tss = int(li[3]) if li[6] == '+' else int(li[4])
        wstart = max(tss - flank, 0)
        for b in range(0, 2 * flank, binsize):
            out.write(f'{chrom}\t{wstart + b}\t{wstart + b + binsize}\t{gid}\n')
