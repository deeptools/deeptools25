# Adjust these if you want
ORGANISM = "human"
PROTOCOL = "chip"
FULL_GTF = False
BINSIZE = 10
UPSTREAM = 500
DOWNSTREAM = 1500

Ntimes = 1
Nthreads = 4


# Do not edit any further
if FULL_GTF:
    GTF = { "human": "regions/homo.v91.full.gtf", "wheat": "regions/triticum.v60.full.gtf" }
else:
    # These were generated via: grep 'transcript_id' full.gtf | shuf | head -n 1000 > sample.gtf
    GTF = { "human": "regions/homo.v91.sample.gtf", "wheat": "regions/triticum.v60.sample.gtf" }

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
        expand("output/benchmark_bamCoverage_bs" + str(BINSIZE) + "_plot_{type}.png", type=["time", "mem"]),
        expand("output/benchmark_bamCompare_bs" + str(BINSIZE) + "_plot_{type}.png",type=["time", "mem"]),
        expand("output/benchmark_computeMatrix_bs" + str(BINSIZE) + "_plot_{type}.png",type=["time", "mem"])


rule bamCoverage2:
    input:
        bam = FILES[ORGANISM + "_" + PROTOCOL]
    output:
        bed = "output/new2.bg"
    benchmark:
        repeat("output/benchmark_bamCoverage2.txt", Ntimes)
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
        bam = FILES[ORGANISM + "_" + PROTOCOL]
    output:
        bed = "output/new1.bg"
    benchmark:
        repeat("output/benchmark_bamCoverage1.txt", Ntimes)
    params:
        binsize = BINSIZE
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        bamCoverage -b {input.bam} -o {output.bed} -of bedgraph -bs {params.binsize} -p {threads}
        """


rule bamCompare2:
    input:
        bam = FILES[ORGANISM + "_" + PROTOCOL]
    output:
        bw = "output/bamCompare2.bw"
    benchmark:
        repeat("output/benchmark_bamCompare2.txt", Ntimes)
    params:
        binsize = BINSIZE
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        bamCompare -b1 {input.bam} -b2 {input.bam} -o {output.bw} -bs {params.binsize} -p {threads} --exactScaling
        """

rule bamCompare1:
    input:
        bam = FILES[ORGANISM + "_" + PROTOCOL],
    output:
        bw = "output/bamCompare1.bw"
    benchmark:
        repeat("output/benchmark_bamCompare1.txt", Ntimes)
    params:
        binsize = BINSIZE
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        bamCompare -b1 {input.bam} -b2 {input.bam} -o {output.bw} -bs {params.binsize} -p {threads} --exactScaling
        """


rule computeMatrix2:
    input:
        bw1 = "output/bamCompare1.bw",
        bw2 = "output/bamCompare2.bw",
        bed = GTF[ORGANISM]
    output:
        npz = "output/test_new2.npz"
    benchmark:
        repeat("output/benchmark_computeMatrix2.txt", Ntimes)
    params:
        binsize = BINSIZE,
        upstream = UPSTREAM,
        downstream = DOWNSTREAM
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        computeMatrix reference-point -S {input.bw1} {input.bw2} -R {input.bed} {input.bed} -o {output.npz} \
            -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero
        """

rule computeMatrix1:
    input:
        bw1 = "output/bamCompare1.bw",
        bw2 = "output/bamCompare2.bw",
        bed = GTF[ORGANISM]
    output:
        npz = "output/test_new1.npz"
    benchmark:
        repeat("output/benchmark_computeMatrix1.txt", Ntimes)
    params:
        binsize = BINSIZE,
        upstream = UPSTREAM,
        downstream = DOWNSTREAM
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        computeMatrix reference-point -S {input.bw1} {input.bw2} -R {input.bed} {input.bed} -o {output.npz} \
            -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero
        """


rule multiBamSummary2:
    input:
        bam = FILES[ORGANISM + "_" + PROTOCOL]
    output:
        npz = "output/mb_summary2.npz",
        outraw = "output/mb_summary2.outraw.tab"
    benchmark:
        repeat("output/benchmark_multiBamSummary2.txt", Ntimes)
    params:
        binsize = BINSIZE * 10000,
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam} {input.bam} -o {output.npz} \
            --outRawCounts {output.outraw} \
            -bs {params.binsize} -p {threads} > /dev/null
        touch {output.npz}
        """

rule multiBamSummary1:
    input:
        bam = FILES[ORGANISM + "_" + PROTOCOL]
    output:
        npz = "output/mb_summary1.npz",
        outraw = "output/mb_summary1.outraw.tab"
    benchmark:
        repeat("output/benchmark_multiBamSummary1.txt", Ntimes)
    params:
        binsize = BINSIZE * 10000,
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam} {input.bam} -o {output.npz} \
            --outRawCounts {output.outraw} \
            -bs {params.binsize} -p {threads} > /dev/null
        touch {output.npz}
        """


rule plot_all_benchmarks:
    input:
        benchmark_bamCoverage1 = "output/benchmark_bamCoverage1.txt",
        benchmark_bamCoverage2 = "output/benchmark_bamCoverage2.txt",
        benchmark_bamCompare1 = "output/benchmark_bamCompare1.txt",
        benchmark_bamCompare2 = "output/benchmark_bamCompare2.txt",
        benchmark_computeMatrix1 = "output/benchmark_computeMatrix1.txt",
        benchmark_computeMatrix2 = "output/benchmark_computeMatrix2.txt",
        benchmark_multiBamSummary1 = "output/benchmark_multiBamSummary1.txt",
        benchmark_multiBamSummary2 = "output/benchmark_multiBamSummary2.txt"
    output:
        bamCoverage_time_plot = "output/benchmark_bamCoverage_bs" + str(BINSIZE) + "_plot_time.png",
        bamCoverage_mem_plot = "output/benchmark_bamCoverage_bs" + str(BINSIZE) + "_plot_mem.png",
        bamCompare_time_plot = "output/benchmark_bamCompare_bs" + str(BINSIZE) + "_plot_time.png",
        bamCompare_mem_plot = "output/benchmark_bamCompare_bs" + str(BINSIZE) + "_plot_mem.png",
        computeMatrix_time_plot = "output/benchmark_computeMatrix_bs" + str(BINSIZE) + "_plot_time.png",
        computeMatrix_mem_plot = "output/benchmark_computeMatrix_bs" + str(BINSIZE) + "_plot_mem.png"
    params:
        bamCoverage_template = "output/benchmark_bamCoverage_bs" + str(BINSIZE) + "_plot.png",
        bamCompare_template = "output/benchmark_bamCompare_bs" + str(BINSIZE) + "_plot.png",
        computeMatrix_template = "output/benchmark_computeMatrix_bs" + str(BINSIZE) + "_plot.png"
    conda: "v4.env.yaml"
    shell:
        """
            plot.py {params.bamCoverage_template} {input.benchmark_bamCoverage1} {input.benchmark_bamCoverage2}
            plot.py {params.bamCompare_template} {input.benchmark_bamCompare1} {input.benchmark_bamCompare2}
            plot.py {params.computeMatrix_template} {input.benchmark_computeMatrix1} {input.benchmark_computeMatrix2}
        """

