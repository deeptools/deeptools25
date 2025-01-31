# Adjust these if you want
ORGANISM = "human"
FULL_GTF = False
BINSIZE = 10
UPSTREAM = 500
DOWNSTREAM = 1500

Ntimes = 1
Nthreads = 12


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
        expand("output/bamCoverage_{organism}_bs{binsize}_plot_{type}.png", 
               organism=ORGANISM, binsize=BINSIZE, type=["time", "mem"]),
        expand("output/bamCompare_{organism}_bs{binsize}_plot_{type}.png", 
               organism=ORGANISM, binsize=BINSIZE, type=["time", "mem"])
        # expand("output/new1_{protocol}.bg", protocol=PROTOCOLS),
        # expand("output/new2_{protocol}.bg", protocol=PROTOCOLS),
        # expand("output/computeMatrix_bs" + str(BINSIZE) + "_plot_{type}.png", type=["time", "mem"])


rule bamCoverage2:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = "output/new2_{protocol}.bg"
    benchmark:
        repeat("output/bamCoverage2_{protocol}.txt", Ntimes)
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
        repeat("output/bamCoverage1_{protocol}.txt", Ntimes)
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
        repeat("output/bamCompare2.txt", Ntimes)
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
        repeat("output/bamCompare1.txt", Ntimes)
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
        repeat("output/computeMatrix2.txt", Ntimes)
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
        repeat("output/computeMatrix1.txt", Ntimes)
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
        repeat("output/multiBamSummary2.txt", Ntimes)
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
        repeat("output/multiBamSummary1.txt", Ntimes)
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
        bamCoverage1 = expand("output/bamCoverage1_{protocol}.txt", protocol=["chip", "rna", "wgs"]),
        bamCoverage2 = expand("output/bamCoverage2_{protocol}.txt", protocol=["chip", "rna", "wgs"]),
        bamCompare1 = "output/bamCompare1.txt",
        bamCompare2 = "output/bamCompare2.txt",
        # computeMatrix1 = "output/computeMatrix1.txt",
        # computeMatrix2 = "output/computeMatrix2.txt"
    output:
        bamCoverage_time_plot = "output/bamCoverage_{ORGANISM}_bs{BINSIZE}_plot_time.png",
        bamCoverage_mem_plot = "output/bamCoverage_{ORGANISM}_bs{BINSIZE}_plot_mem.png",
        bamCompare_time_plot = "output/bamCompare_{ORGANISM}_bs{BINSIZE}_plot_time.png",
        bamCompare_mem_plot = "output/bamCompare_{ORGANISM}_bs{BINSIZE}_plot_mem.png",
        # computeMatrix_time_plot = "output/computeMatrix_bs{BINSIZE}_plot_time.png",
        # computeMatrix_mem_plot = "output/computeMatrix_bs{BINSIZE}_plot_mem.png"
    params:
        bamCoverage_template = "output/bamCoverage_{ORGANISM}_bs{BINSIZE}_plot.png",
        bamCompare_template = "output/bamCompare_{ORGANISM}_bs{BINSIZE}_plot.png",
        # computeMatrix_template = "output/computeMatrix_bs{BINSIZE}_plot.png"
    shell:
        """
        plot.py {params.bamCoverage_template} \
            {input.benchmark_bamCoverage1[0]},{input.benchmark_bamCoverage1[1]},{input.benchmark_bamCoverage1[2]} \
            {input.benchmark_bamCoverage2[0]},{input.benchmark_bamCoverage2[1]},{input.benchmark_bamCoverage2[2]}
        plot.py {params.bamCompare_template} {input.benchmark_bamCompare1} {input.benchmark_bamCompare2}
        # plot.py {params.computeMatrix_template} {input.computeMatrix1} {input.computeMatrix2}
        """