import subprocess as sp
from multiprocessing import Pool
from pathlib import Path
import shutil

its = snakemake.params.its


def fqdump(tup):
    sample, sra, ofq = tup
    ofp = Path(ofq)
    ofp.mkdir(exist_ok=True)
    if not list( ofp.glob(f"{sample}*.fastq.gz") ):
        if not Path(sra, sra + '.sra').exists():
          prefc = ['prefetch', '-X', '500G', sra]
          sp.run(prefc)
        if not list( ofp.glob(f"{sra}*.fastq.gz") ):
            splitc = ['fastq-dump', '--split-3', '-O', str(ofp), f"{sra}/{sra}.sra"]
            sp.run(splitc)
            # R1 / R2
            (ofp/f"{sra}_1.fastq").rename( (ofp/f"{sra}_R1.fastq") )
            (ofp/f"{sra}_2.fastq").rename( (ofp/f"{sra}_R2.fastq") )
            compressc = ['pigz', '-9', str(ofp/f"{sra}_R1.fastq")]
            sp.run(compressc)
            compressc = ['pigz', '-9', str(ofp/f"{sra}_R2.fastq")]
            sp.run(compressc)
        (ofp/f"{sra}_R1.fastq.gz").rename( (ofp/f"{sample}_R1.fastq.gz") )
        (ofp/f"{sra}_R2.fastq.gz").rename( (ofp/f"{sample}_R2.fastq.gz") )
    if Path(sra).exists() and Path(sra).is_dir():
        shutil.rmtree(sra)

with Pool(snakemake.threads) as p:
    p.map(fqdump, its)
