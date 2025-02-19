# Adjust these if you want
ORGANISM = "homo"
GTF = "regions/homo.v91.sample25k.gtf"

BinSizes = {
    "bamCoverage": 10,
    "bamCompare": 100,
    "computeMatrix": 300,
    "multiBamSummary": 5000,
}

Ntimes = 3
Nthreads = 10


# Do not edit any further
import platform
system = platform.system()
if system == "Linux":
    timeCmd = "/usr/bin/time -v"
elif system == "Darwin":
    timeCmd = "/usr/bin/time -al"
else:
    raise ValueError(f"Unsupported platform: {system}")

PROTOCOLS = ["chip", "rna", "wgs"]
FILES = {
    "homo_chip": "zenodo/human_chip_SRR28592124.bam",
    "homo_rna": "zenodo/human_rna_SRR28012902.bam",
    "homo_wgs": "zenodo/human_wgs_SRR15494527.bam",
    "triticum_chip": "zenodo/triticum_chip_SRR1686799.bam",
    "triticum_rna": "zenodo/triticum_rna_SRR27822150.bam",
    "triticum_wgs": "zenodo/triticum_wgs_SRR27887047.bam"
}

# For human we have a more intensive computeMatrix comparison
# For Triticum, we'll simply repeat the same BW a couple of times. For multiBamSummary (both sp.) we'll do something similar (repeat same input many times.)
humanBigwigs = "zenodo/bigwigs/human_chip_SRR28592124.bw zenodo/bigwigs/human_chip_SRR28592125.bw zenodo/bigwigs/human_chip_SRR28592131.bw zenodo/bigwigs/human_chip_SRR28592132.bw zenodo/bigwigs/human_rna_SRR28012902.bw zenodo/bigwigs/human_rna_SRR28012903.bw zenodo/bigwigs/human_rna_SRR28012904.bw zenodo/bigwigs/human_rna_SRR28012905.bw zenodo/bigwigs/human_wgs_SRR15494527.bw"

# TODO: if we were to use triticumBigwigs, adjut README.md accordingly,
# ├── triticum_chip_SRR1686799_fw.bw
# ├── triticum_chip_SRR1686799_mapq10.bw
# ├── triticum_chip_SRR1686799_markdup.bw
# ├── triticum_chip_SRR1686799_nodup.bw
# ├── triticum_chip_SRR1686799_reverse.bw

rule all:
    input:
        expand("output/bamCoverage_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["bamCoverage"], type=["time", "mem"]),
        expand("output/bamCompare_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["bamCompare"], type=["time", "mem"]),
        expand("output/computeMatrix_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["computeMatrix"], type=["time", "mem"]),
        expand("output/multiBamSummary_{organism}_bs{binsize}_{type}.png",
               organism=ORGANISM, binsize=BinSizes["multiBamSummary"], type=["time", "mem"])


rule bamCoverage2:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = "output/bamCoverage2_{protocol}.bg",
        iter_file = "output/benchmark_iteration_bamCoverage2_{protocol}.txt"
    benchmark:
        repeat(f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes["bamCoverage"]}_{{protocol}}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCoverage"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        {timeCmd} bamCoverage -b {input.bam} -o {output.bed} -of bedgraph \
            -bs {params.binsize} -p {threads} \
                >logs/bamCoverage2_{wildcards.protocol}_${{curr_iter}}.txt 2>&1
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule bamCoverage1:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = "output/bamCoverage1_{protocol}.bg",
        iter_file = "output/benchmark_iteration_bamCoverage1_{protocol}.txt"
    benchmark:
        repeat(f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes["bamCoverage"]}_{{protocol}}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCoverage"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        {timeCmd} bamCoverage -b {input.bam} -o {output.bed} -of bedgraph \
            -bs {params.binsize} -p {threads} \
                >logs/bamCoverage1_{wildcards.protocol}_${{curr_iter}}.txt 2>&1
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule bamCompare2:
    input:
        bam1 = lambda wildcards: FILES[ORGANISM + "_chip"],
        bam2 = lambda wildcards: FILES[ORGANISM + "_wgs"]
    output:
        bw = "output/bamCompare2.bw",
        iter_file = "output/benchmark_iteration_bamCompare2.txt"
    benchmark:
        repeat(f"logs/bamCompare2_{ORGANISM}_bs{BinSizes["bamCompare"]}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCompare"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o {output.bw} -bs {params.binsize} -p {threads} \
                >logs/bamCompare2_${{curr_iter}}.txt 2>&1
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule bamCompare1:
    input:
        bam1 = lambda wildcards: FILES[ORGANISM + "_chip"],
        bam2 = lambda wildcards: FILES[ORGANISM + "_wgs"]
    output:
        bw = "output/bamCompare1.bw",
        iter_file = "output/benchmark_iteration_bamCompare1.txt"
    benchmark:
        repeat(f"logs/bamCompare1_{ORGANISM}_bs{BinSizes["bamCompare"]}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCompare"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o {output.bw} -bs {params.binsize} -p {threads} \
                >logs/bamCompare1_${{curr_iter}}.txt 2>&1
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule computeMatrix2:
    input:
        bw2 = "output/bamCompare2.bw",
        gtf = GTF
    output:
        npz = "output/computeMatrix2.npz",
        iter_file = "output/benchmark_iteration_computeMatrix2.txt"
    benchmark:
        repeat(f"logs/computeMatrix2_{ORGANISM}_bs{BinSizes["computeMatrix"]}.txt", Ntimes)
    params:
        binsize = BinSizes["computeMatrix"],
        upstream = 2 * BinSizes["computeMatrix"],
        downstream = 2 * BinSizes["computeMatrix"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        if [ "{ORGANISM}" = "homo" ]; then
          {timeCmd} computeMatrix reference-point \
              -S {input.bw2} {humanBigwigs} \
              -R {input.gtf} -o {output.npz} -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  >logs/computeMatrix2_${{curr_iter}}.txt 2>&1
        else
          {timeCmd} computeMatrix reference-point \
              -S {input.bw2} {input.bw2} {input.bw2} {input.bw2} {input.bw2} {input.bw2} {input.bw2} {input.bw2} {input.bw2} {input.bw2} \
              -R {input.gtf} -o {output.npz} -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  >logs/computeMatrix2_${{curr_iter}}.txt 2>&1
        fi
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule computeMatrix1:
    input:
        bw1 = "output/bamCompare1.bw",
        gtf = GTF
    output:
        npz = "output/computeMatrix1.npz",
        iter_file = "output/benchmark_iteration_computeMatrix1.txt"
    benchmark:
        repeat(f"logs/computeMatrix1_{ORGANISM}_bs{BinSizes["computeMatrix"]}.txt", Ntimes)
    params:
        binsize = BinSizes["computeMatrix"],
        upstream = 2 * BinSizes["computeMatrix"],
        downstream = 2 * BinSizes["computeMatrix"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        if [ "{ORGANISM}" = "homo" ]; then
          {timeCmd} computeMatrix reference-point \
              -S {input.bw1} {humanBigwigs} \
              -R {input.gtf} -o {output.npz} -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  >logs/computeMatrix1_${{curr_iter}}.txt 2>&1
        else
          {timeCmd} computeMatrix reference-point \
              -S {input.bw1} {input.bw1} {input.bw1} {input.bw1} {input.bw1} {input.bw1} {input.bw1} {input.bw1} {input.bw1} {input.bw1} \
              -R {input.gtf} -o {output.npz} -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  >logs/computeMatrix1_${{curr_iter}}.txt 2>&1
        fi
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule multiBamSummary2:
    input:
        bam = FILES[ORGANISM + "_" + "wgs"]
    output:
        npz = "output/multiBamSummary2.npz",
        outraw = "output/multiBamSummary2.outraw.tab",
        iter_file = "output/benchmark_iteration_multiBamSummary2.txt"
    benchmark:
        repeat(f"logs/multiBamSummary2_{ORGANISM}_bs{BinSizes["multiBamSummary"]}.txt", Ntimes)
    params:
        binsize = BinSizes["multiBamSummary"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        if [ "{ORGANISM}" = "homo" ]; then
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} {input.bam} {input.bam} \
                -o {output.npz} --outRawCounts {output.outraw} -bs {params.binsize} -p {threads} \
                    >logs/multiBamSummary2_${{curr_iter}}.txt 2>&1
        else
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} \
                -o {output.npz} --outRawCounts {output.outraw} -bs {params.binsize} -p {threads} \
                    >logs/multiBamSummary2_${{curr_iter}}.txt 2>&1
        fi
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule multiBamSummary1:
    input:
        bam = FILES[ORGANISM + "_" + "wgs"]
    output:
        npz = "output/multiBamSummary1.npz",
        outraw = "output/multiBamSummary1.outraw.tab",
        iter_file = "output/benchmark_iteration_multiBamSummary1.txt"
    benchmark:
        repeat(f"logs/multiBamSummary1_{ORGANISM}_bs{BinSizes["multiBamSummary"]}.txt", Ntimes)
    params:
        binsize = BinSizes["multiBamSummary"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        if [ "{ORGANISM}" = "homo" ]; then
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} {input.bam} {input.bam} \
                -o {output.npz} --outRawCounts {output.outraw} -bs {params.binsize} -p {threads} \
                    >logs/multiBamSummary1_${{curr_iter}}.txt 2>&1
        else
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} \
                -o {output.npz} --outRawCounts {output.outraw} -bs {params.binsize} -p {threads} \
                    >logs/multiBamSummary1_${{curr_iter}}.txt 2>&1
        fi
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """


rule process_results:
    input:
        bamCoverage1_chip = f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes["bamCoverage"]}_chip.txt",
        bamCoverage2_chip = f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes["bamCoverage"]}_chip.txt",
        bamCoverage1_rna = f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes["bamCoverage"]}_rna.txt",
        bamCoverage2_rna = f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes["bamCoverage"]}_rna.txt",
        bamCoverage1_wgs = f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes["bamCoverage"]}_wgs.txt",
        bamCoverage2_wgs = f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes["bamCoverage"]}_wgs.txt",
        bamCompare1 = f"logs/bamCompare1_{ORGANISM}_bs{BinSizes["bamCompare"]}.txt",
        bamCompare2 = f"logs/bamCompare2_{ORGANISM}_bs{BinSizes["bamCompare"]}.txt",
        computeMatrix1 = f"logs/computeMatrix1_{ORGANISM}_bs{BinSizes["computeMatrix"]}.txt",
        computeMatrix2 = f"logs/computeMatrix2_{ORGANISM}_bs{BinSizes["computeMatrix"]}.txt",
        multiBamSummary1 = f"logs/multiBamSummary1_{ORGANISM}_bs{BinSizes["multiBamSummary"]}.txt",
        multiBamSummary2 = f"logs/multiBamSummary2_{ORGANISM}_bs{BinSizes["multiBamSummary"]}.txt"
    output:
        bamCoverage_time_plot = f"output/bamCoverage_{ORGANISM}_bs{BinSizes["bamCoverage"]}_time.png",
        bamCoverage_mem_plot = f"output/bamCoverage_{ORGANISM}_bs{BinSizes["bamCoverage"]}_mem.png",
        bamCompare_time_plot = f"output/bamCompare_{ORGANISM}_bs{BinSizes["bamCompare"]}_time.png",
        bamCompare_mem_plot = f"output/bamCompare_{ORGANISM}_bs{BinSizes["bamCompare"]}_mem.png",
        computeMatrix_time_plot = f"output/computeMatrix_{ORGANISM}_bs{BinSizes["computeMatrix"]}_time.png",
        computeMatrix_mem_plot = f"output/computeMatrix_{ORGANISM}_bs{BinSizes["computeMatrix"]}_mem.png",
        multiBamSummary_time_plot = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes["multiBamSummary"]}_time.png",
        multiBamSummary_mem_plot = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes["multiBamSummary"]}_mem.png"
    shell:
        """
        python3 present_results.py output/bamCoverage_{ORGANISM}_bs{BinSizes["bamCoverage"]}.png \
            {input.bamCoverage1_chip},{input.bamCoverage1_rna},{input.bamCoverage1_wgs} \
            {input.bamCoverage2_chip},{input.bamCoverage2_rna},{input.bamCoverage2_wgs}
        python3 present_results.py output/bamCompare_{ORGANISM}_bs{BinSizes["bamCompare"]}.png {input.bamCompare1} {input.bamCompare2}
        python3 present_results.py output/computeMatrix_{ORGANISM}_bs{BinSizes["computeMatrix"]}.png {input.computeMatrix1} {input.computeMatrix2}
        python3 present_results.py output/multiBamSummary_{ORGANISM}_bs{BinSizes["multiBamSummary"]}.png {input.multiBamSummary1} {input.multiBamSummary2}
        """
