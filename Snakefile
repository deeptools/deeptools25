# Adjust these if you want
ORGANISM = "human"
FULL_GTF = False
BINSIZE = 10
UPSTREAM = 500
DOWNSTREAM = 1500

Ntimes = 1
Nthreads = 16


# Do not edit any further
if FULL_GTF:
    GTF = { "human": "regions/homo.v91.full.gtf", "wheat": "regions/triticum.v60.full.gtf" }
else:
    GTF = { "human": "regions/homo.v91.sample.gtf", "wheat": "regions/triticum.v60.sample.gtf" }

PROTOCOLS = ["chip", "rna", "wgs"]
FILES = {
    "human_chip": "zenodo/human_chip_SRR28592124.bam",
    "human_rna": "zenodo/human_rna_SRR28012902.bam",
    "human_wgs": "zenodo/human_wgs_SRR15494527.bam",
    "triticum_chip": "zenodo/triticum_chip_SRR1686799.bam",
    "triticum_rna": "zenodo/triticum_rna_SRR27822150.bam",
    "triticum_wgs": "zenodo/triticum_wgs_SRR27887047.bam"
}


rule all:
    input:
        expand("output/bamCoverage_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BINSIZE, type=["time", "mem"]),
        expand("output/bamCompare_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BINSIZE, type=["time", "mem"]),
        expand("output/computeMatrix_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BINSIZE, type=["time", "mem"]),
        expand("output/multiBamSummary_{organism}_bs{binsize}_{type}.png",
               organism=ORGANISM, binsize=BINSIZE, type=["time", "mem"])

rule bamCoverage2:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = "output/new2_{protocol}.bg"
    benchmark:
        repeat(f"output/bamCoverage2_{ORGANISM}_bs{BINSIZE}_{wildcards.protocol}.txt", Ntimes)
    params:
        binsize = BINSIZE
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        bamCoverage -b {input.bam} -o {output.bed} -of bedgraph -bs {params.binsize} -p {threads}
        """

rule bamCoverage1:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = "output/new1_{protocol}.bg"
    benchmark:
        repeat(f"output/bamCoverage1_{ORGANISM}_bs{BINSIZE}_{wildcards.protocol}.txt", Ntimes)
    params:
        binsize = BINSIZE
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        bamCoverage -b {input.bam} -o {output.bed} -of bedgraph -bs {params.binsize} -p {threads} > /dev/null
        """



rule bamCompare2:
    input:
        bam1 = FILES[ORGANISM + "_chip"],
        bam2 = FILES[ORGANISM + "_wgs"]
    output:
        bw = "output/bamCompare2.bw"
    benchmark:
        repeat(f"output/bamCompare2_{ORGANISM}_bs{BINSIZE}.txt", Ntimes)
    params:
        binsize = BINSIZE
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        bamCompare -b1 {input.bam1} -b2 {input.bam2} -o {output.bw} -bs {params.binsize} -p {threads} --exactScaling
        """

rule bamCompare1:
    input:
        bam1 = FILES[ORGANISM + "_chip"],
        bam2 = FILES[ORGANISM + "_wgs"]
    output:
        bw = "output/bamCompare1.bw"
    benchmark:
        repeat(f"output/bamCompare1_{ORGANISM}_bs{BINSIZE}.txt", Ntimes)
    params:
        binsize = BINSIZE
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        bamCompare -b1 {input.bam1} -b2 {input.bam2} -o {output.bw} -bs {params.binsize} -p {threads} --exactScaling > /dev/null
        """


rule computeMatrix2:
    input:
        bw2 = "output/bamCompare2.bw",
        bed = GTF[ORGANISM]
    output:
        npz = "output/test_new2.npz"
    benchmark:
        repeat(f"output/computeMatrix2_{ORGANISM}_bs{BINSIZE}.txt", Ntimes)
    params:
        binsize = BINSIZE,
        upstream = UPSTREAM,
        downstream = DOWNSTREAM
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        computeMatrix reference-point -S {input.bw2} {input.bw2} -R {input.bed} {input.bed} -o {output.npz} \
            -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero
        """

rule computeMatrix1:
    input:
        bw1 = "output/bamCompare1.bw",
        bed = GTF[ORGANISM]
    output:
        npz = "output/test_new1.npz"
    benchmark:
        repeat(f"output/computeMatrix1_{ORGANISM}_bs{BINSIZE}.txt", Ntimes)
    params:
        binsize = BINSIZE,
        upstream = UPSTREAM,
        downstream = DOWNSTREAM
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        computeMatrix reference-point -S {input.bw1} {input.bw1} -R {input.bed} {input.bed} -o {output.npz} \
            -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero > /dev/null
        """


rule multiBamSummary2:
    input:
        bam = FILES[ORGANISM + "_" + "wgs"]
    output:
        npz = "output/mb_summary2.npz",
        outraw = "output/mb_summary2.outraw.tab"
    benchmark:
        repeat(f"output/multiBamSummary2_{ORGANISM}_bs{BINSIZE}.txt", Ntimes)
    params:
        binsize = BINSIZE * 10000,
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam} {input.bam} -o {output.npz} \
            --outRawCounts {output.outraw} \
            -bs {params.binsize} -p {threads}
        """

rule multiBamSummary1:
    input:
        bam = FILES[ORGANISM + "_" + "wgs"]
    output:
        npz = "output/mb_summary1.npz",
        outraw = "output/mb_summary1.outraw.tab"
    benchmark:
        repeat(f"output/multiBamSummary1_{ORGANISM}_bs{BINSIZE}.txt", Ntimes)
    params:
        binsize = BINSIZE * 10000,
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam} {input.bam} -o {output.npz} \
            --outRawCounts {output.outraw} \
            -bs {params.binsize} -p {threads} > /dev/null
        """


rule plot_all_benchmarks:
    input:
        bamCoverage1_chip = f"output/bamCoverage1_{ORGANISM}_bs{BINSIZE}_chip.txt",
        bamCoverage1_rna = f"output/bamCoverage1_{ORGANISM}_bs{BINSIZE}_rna.txt",
        bamCoverage1_wgs = f"output/bamCoverage1_{ORGANISM}_bs{BINSIZE}_wgs.txt",
        bamCoverage2_chip = f"output/bamCoverage2_{ORGANISM}_bs{BINSIZE}_chip.txt",
        bamCoverage2_rna = f"output/bamCoverage2_{ORGANISM}_bs{BINSIZE}_rna.txt",
        bamCoverage2_wgs = f"output/bamCoverage2_{ORGANISM}_bs{BINSIZE}_wgs.txt",
        bamCompare1 = f"output/bamCompare1_{ORGANISM}_bs{BINSIZE}.txt",
        bamCompare2 = f"output/bamCompare2_{ORGANISM}_bs{BINSIZE}.txt",
        computeMatrix1 = f"output/computeMatrix1_{ORGANISM}_bs{BINSIZE}.txt",
        computeMatrix2 = f"output/computeMatrix2_{ORGANISM}_bs{BINSIZE}.txt",
        multiBamSummary1 = f"output/multiBamSummary1_{ORGANISM}_bs{BINSIZE}.txt",
        multiBamSummary2 = f"output/multiBamSummary2_{ORGANISM}_bs{BINSIZE}.txt"
    output:
        bamCoverage_time_plot = f"output/bamCoverage_{ORGANISM}_bs{BINSIZE}_time.png",
        bamCoverage_mem_plot = f"output/bamCoverage_{ORGANISM}_bs{BINSIZE}_mem.png",
        bamCompare_time_plot = f"output/bamCompare_{ORGANISM}_bs{BINSIZE}_time.png",
        bamCompare_mem_plot = f"output/bamCompare_{ORGANISM}_bs{BINSIZE}_mem.png",
        computeMatrix_time_plot = f"output/computeMatrix_{ORGANISM}_bs{BINSIZE}_time.png",
        computeMatrix_mem_plot = f"output/computeMatrix_{ORGANISM}_bs{BINSIZE}_mem.png",
        multiBamSummary_time_plot = f"output/multiBamSummary_{ORGANISM}_bs{BINSIZE}_time.png",
        multiBamSummary_mem_plot = f"output/multiBamSummary_{ORGANISM}_bs{BINSIZE}_mem.png"
    shell:
        """
        python3 plot.py output/bamCoverage_{ORGANISM}_bs{BINSIZE}.png \
            {input.bamCoverage1_chip},{input.bamCoverage1_rna},{input.bamCoverage1_wgs} \
            {input.bamCoverage2_chip},{input.bamCoverage2_rna},{input.bamCoverage2_wgs}
        python3 plot.py output/bamCompare_{ORGANISM}_bs{BINSIZE}.png {input.bamCompare1} {input.bamCompare2}
        python3 plot.py output/computeMatrix_{ORGANISM}_bs{BINSIZE}.png {input.computeMatrix1} {input.computeMatrix2}
        python3 plot.py output/multiBamSummary_{ORGANISM}_bs{BINSIZE}.png {input.multiBamSummary1} {input.multiBamSummary2}
        """
