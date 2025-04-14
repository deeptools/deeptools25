# Parameters for the workflow, with defaults to be overridden with snakemake, e.g. `--config organism=triticum`
Ntimes = config.get("ntimes", 3)
Nthreads = config.get("nthreads", 16)
ORGANISM = config.get("organism", "homo")

# These are customizable too
BinSizes = {
    "bamCoverage": 10,
    "bamCompare": 100,
    "computeMatrix": 300,
    "multiBamSummary": 500,
}

# Helper fn. to keep logs of failed jobs
shell.prefix("""
function on_error() {{ 
    cp $1 $1.failed.$(date +%Y%m%d_%H%M%S)
    echo "ERROR: Command failed, log saved as $1.failed.$(date +%Y%m%d_%H%M%S)" >&2
    
    # If we have an iteration file, decrement it so the job will be retried
    if [ -n "$2" ] && [ -f "$2" ]; then
        curr_iter=$(cat $2)
        if [ $curr_iter -gt 1 ]; then
            new_iter=$((curr_iter - 1))
            echo $new_iter > $2
            echo "Decremented iteration counter in $2 from $curr_iter to $new_iter" >&2
        fi
    fi
    
    return 1  # This will propagate the error
}}; 
"""
)

# Set the time command based on the platform
import platform
system = platform.system()
if system == "Linux":
    timeCmd = "/usr/bin/time -v"
elif system == "Darwin":
    timeCmd = "/usr/bin/time -al"
else:
    raise ValueError(f"Unsupported platform: {system}")

# Input files
PROTOCOLS = ["chip", "rna", "wgs"]
FILES = {
    "homo_chip": "zenodo/human_chip_SRR28592124.bam",
    "homo_rna": "zenodo/human_rna_SRR28012902.bam",
    "homo_wgs": "zenodo/human_wgs_SRR15494527.bam",
    "triticum_chip": "zenodo/triticum_chip_SRR1686799.bam",
    "triticum_rna": "zenodo/triticum_rna_SRR27822150.bam",
    "triticum_wgs": "zenodo/triticum_wgs_SRR27887047.bam"
}

# We have a more intensive computeMatrix comparison
humanBigwigs = "zenodo/bigwigs/human_chip_SRR28592124.bw zenodo/bigwigs/human_chip_SRR28592125.bw zenodo/bigwigs/human_chip_SRR28592131.bw zenodo/bigwigs/human_chip_SRR28592132.bw zenodo/bigwigs/human_rna_SRR28012902.bw zenodo/bigwigs/human_rna_SRR28012903.bw zenodo/bigwigs/human_rna_SRR28012904.bw zenodo/bigwigs/human_rna_SRR28012905.bw zenodo/bigwigs/human_wgs_SRR15494527.bw"
wheatBigwigs = " zenodo/bigwigs/triticum_chip_SRR1686799_mapq10.bw zenodo/bigwigs/triticum_chip_SRR1686799_nodup.bw zenodo/bigwigs/triticum_rna_SRR27822150_mapq10.bw zenodo/bigwigs/triticum_wgs_mapq20.bw"
# For multiBamSummary (both sp.) we'll do something similar (repeat same input many times.)

# The number of transcripts is hardcoded in filenames here, sorry.
if ORGANISM == "homo":
    GTF = f"regions/{ORGANISM}.v91.full.gtf"
elif ORGANISM == "triticum":
    GTF = f"regions/{ORGANISM}.v60.full.gtf"
else:
    raise ValueError(f"Unsupported organism: {ORGANISM}")

rule all:
    input: "output/report.html"

#rule downsample_gtf:
#    input: GTF.replace('sample25k', 'full')
#    output: GTF
#    conda: "extras.env.yaml"
#    log: "logs/downsample_gtf.log"
#    params:
#        transcript_count = 25000
#    shell:
#        """
#        # Check if input file exists and has sufficient entries
#        if [ ! -f {input} ]; then
#            echo "Error: Input GTF file {input} does not exist" > {log}
#            exit 1
#        fi
#        transcript_count=$(grep -c 'transcript_id' {input})
#        if [ $transcript_count -lt {params.transcript_count} ]; then
#            echo "Warning: Input GTF doesn't have enough transcripts ($transcript_count)" > {log}
#        fi
#        
#        # Perform downsampling
#        grep 'transcript_id' {input} | shuf | head -n {params.transcript_count} | bedtools sort -i - > {output}
#        """

rule bamCoverage2:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = directory("output/bamCoverage2_{protocol}"),
        iter_file = "output/iter_count_bamCoverage2_{protocol}.txt",
        done = expand("output/bamCoverage2_{{protocol}}_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_bamCoverage2_{ORGANISM}_bs{BinSizes['bamCoverage']}_{{protocol}}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCoverage"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        mkdir -p {output.bed}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bed}/iter_${{curr_iter}}.bg"
        log_file="logs/bamCoverage2_{wildcards.protocol}_${{curr_iter}}.txt"
        
        {timeCmd} bamCoverage -b {input.bam} -o $out_file -of bedgraph \
            -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
                
        # Create done marker file
        touch output/bamCoverage2_{wildcards.protocol}_done_${{curr_iter}}.txt
                
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule bamCoverage1:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = directory("output/bamCoverage1_{protocol}"),
        iter_file = "output/iter_count_bamCoverage1_{protocol}.txt",
        done = expand("output/bamCoverage1_{{protocol}}_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_bamCoverage1_{ORGANISM}_bs{BinSizes['bamCoverage']}_{{protocol}}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCoverage"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        mkdir -p {output.bed}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bed}/iter_${{curr_iter}}.bg"
        log_file="logs/bamCoverage1_{wildcards.protocol}_${{curr_iter}}.txt"
        
        {timeCmd} bamCoverage -b {input.bam} -o $out_file -of bedgraph \
            -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
                
        # Create done marker file
        touch output/bamCoverage1_{wildcards.protocol}_done_${{curr_iter}}.txt
                
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule bamCompare2:
    input:
        bam1 = lambda wildcards: FILES[ORGANISM + "_chip"],
        bam2 = lambda wildcards: FILES[ORGANISM + "_wgs"]
    output:
        bw = directory("output/bamCompare2"),
        iter_file = "output/iter_count_bamCompare2.txt",
        done = expand("output/bamCompare2_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_bamCompare2_{ORGANISM}_bs{BinSizes['bamCompare']}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCompare"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        mkdir -p {output.bw}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bw}/iter_${{curr_iter}}.bw"
        log_file="logs/bamCompare2_${{curr_iter}}.txt"
        
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o $out_file -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
                
        # Create done marker file
        touch output/bamCompare2_done_${{curr_iter}}.txt
                
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule bamCompare1:
    input:
        bam1 = lambda wildcards: FILES[ORGANISM + "_chip"],
        bam2 = lambda wildcards: FILES[ORGANISM + "_wgs"]
    output:
        bw = directory("output/bamCompare1"),
        iter_file = "output/iter_count_bamCompare1.txt",
        done = expand("output/bamCompare1_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_bamCompare1_{ORGANISM}_bs{BinSizes['bamCompare']}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCompare"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        mkdir -p {output.bw}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bw}/iter_${{curr_iter}}.bw"
        log_file="logs/bamCompare1_${{curr_iter}}.txt"
        
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o $out_file -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
                
        # Create done marker file
        touch output/bamCompare1_done_${{curr_iter}}.txt
                
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule computeMatrix2:
    input:
        bw_dir = "output/bamCompare2",
        bw_done = expand("output/bamCompare2_done_{iter}.txt", iter=range(1, Ntimes+1)),
        gtf = GTF
    output:
        npz = directory("output/computeMatrix2"),
        iter_file = "output/iter_count_computeMatrix2.txt",
        done = expand("output/computeMatrix2_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_computeMatrix2_{ORGANISM}_bs{BinSizes['computeMatrix']}.txt", Ntimes)
    params:
        binsize = BinSizes["computeMatrix"],
        upstream = 2 * BinSizes["computeMatrix"],
        downstream = 2 * BinSizes["computeMatrix"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        mkdir -p {output.npz}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.npz}/iter_${{curr_iter}}.npz"
        log_file="logs/computeMatrix2_${{curr_iter}}.txt"
        
        # Try to match iterations, otherwise use first file
        if [ -f "{input.bw_dir}/iter_${{curr_iter}}.bw" ]; then
            input_bw="{input.bw_dir}/iter_${{curr_iter}}.bw"
        else
            input_bw=$(ls {input.bw_dir}/iter_*.bw | head -n 1)
        fi
        
        if [ "{ORGANISM}" = "homo" ]; then
          {timeCmd} computeMatrix reference-point --verbose \
              -S $input_bw {humanBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        else
          {timeCmd} computeMatrix reference-point --verbose \
              -S $input_bw {wheatBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        fi
        
        # Create done marker file
        touch output/computeMatrix2_done_${{curr_iter}}.txt
        
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule computeMatrix1:
    input:
        bw_dir = "output/bamCompare1",
        bw_done = expand("output/bamCompare1_done_{iter}.txt", iter=range(1, Ntimes+1)),
        gtf = GTF
    output:
        npz = directory("output/computeMatrix1"),
        iter_file = "output/iter_count_computeMatrix1.txt",
        done = expand("output/computeMatrix1_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_computeMatrix1_{ORGANISM}_bs{BinSizes['computeMatrix']}.txt", Ntimes)
    params:
        binsize = BinSizes["computeMatrix"],
        upstream = 2 * BinSizes["computeMatrix"],
        downstream = 2 * BinSizes["computeMatrix"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        mkdir -p {output.npz}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.npz}/iter_${{curr_iter}}.npz"
        log_file="logs/computeMatrix1_${{curr_iter}}.txt"
        
        # Try to match iterations, otherwise use first file
        if [ -f "{input.bw_dir}/iter_${{curr_iter}}.bw" ]; then
            input_bw="{input.bw_dir}/iter_${{curr_iter}}.bw"
        else
            input_bw=$(ls {input.bw_dir}/iter_*.bw | head -n 1)
        fi
        
        if [ "{ORGANISM}" = "homo" ]; then
          {timeCmd} computeMatrix reference-point \
              -S $input_bw {humanBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        else
          {timeCmd} computeMatrix reference-point \
              -S $input_bw {wheatBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        fi
        
        # Create done marker file
        touch output/computeMatrix1_done_${{curr_iter}}.txt
        
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule multiBamSummary2:
    input:
        bam = FILES[ORGANISM + "_" + "wgs"]
    output:
        npz = directory("output/multiBamSummary2"),
        iter_file = "output/iter_count_multiBamSummary2.txt",
        done = expand("output/multiBamSummary2_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_multiBamSummary2_{ORGANISM}_bs{BinSizes['multiBamSummary']}.txt", Ntimes)
    params:
        binsize = BinSizes["multiBamSummary"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        mkdir -p {output.npz}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_npz="{output.npz}/iter_${{curr_iter}}.npz"
        out_raw="{output.npz}/iter_${{curr_iter}}.outraw.tab"
        log_file="logs/multiBamSummary2_${{curr_iter}}.txt"
        
        if [ "{ORGANISM}" = "homo" ]; then
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        else
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        fi
        
        # Create done marker file
        touch output/multiBamSummary2_done_${{curr_iter}}.txt
        
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule multiBamSummary1:
    input:
        bam = FILES[ORGANISM + "_" + "wgs"]
    output:
        npz = directory("output/multiBamSummary1"),
        iter_file = "output/iter_count_multiBamSummary1.txt",
        done = expand("output/multiBamSummary1_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/metrics_multiBamSummary1_{ORGANISM}_bs{BinSizes['multiBamSummary']}.txt", Ntimes)
    params:
        binsize = BinSizes["multiBamSummary"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        mkdir -p {output.npz}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_npz="{output.npz}/iter_${{curr_iter}}.npz"
        out_raw="{output.npz}/iter_${{curr_iter}}.outraw.tab"
        log_file="logs/multiBamSummary1_${{curr_iter}}.txt"
        
        if [ "{ORGANISM}" = "homo" ]; then
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        else
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error "$log_file" "{output.iter_file}" && exit 1; }}
        fi
        
        # Create done marker file
        touch output/multiBamSummary1_done_${{curr_iter}}.txt
        
        curr_iter=$((curr_iter + 1))
        echo $curr_iter > {output.iter_file}
        """

rule process_results:
    input:
        # Make sure to wait for all the done markers
        bamCoverage1_done = expand("output/bamCoverage1_{protocol}_done_{iter}.txt", 
                                  protocol=PROTOCOLS, iter=range(1, Ntimes+1)),
        bamCoverage2_done = expand("output/bamCoverage2_{protocol}_done_{iter}.txt", 
                                  protocol=PROTOCOLS, iter=range(1, Ntimes+1)),
        bamCompare1_done = expand("output/bamCompare1_done_{iter}.txt", iter=range(1, Ntimes+1)),
        bamCompare2_done = expand("output/bamCompare2_done_{iter}.txt", iter=range(1, Ntimes+1)),
        computeMatrix1_done = expand("output/computeMatrix1_done_{iter}.txt", iter=range(1, Ntimes+1)),
        computeMatrix2_done = expand("output/computeMatrix2_done_{iter}.txt", iter=range(1, Ntimes+1)),
        multiBamSummary1_done = expand("output/multiBamSummary1_done_{iter}.txt", iter=range(1, Ntimes+1)),
        multiBamSummary2_done = expand("output/multiBamSummary2_done_{iter}.txt", iter=range(1, Ntimes+1)),
        
        # Benchmark files (use the paths generated by Snakemake)
        bamCoverage1_chip = expand("logs/metrics_bamCoverage1_{organism}_bs{binsize}_chip.txt", 
                                  organism=ORGANISM, binsize=BinSizes["bamCoverage"]),
        bamCoverage2_chip = expand("logs/metrics_bamCoverage2_{organism}_bs{binsize}_chip.txt", 
                                  organism=ORGANISM, binsize=BinSizes["bamCoverage"]),
        bamCoverage1_rna = expand("logs/metrics_bamCoverage1_{organism}_bs{binsize}_rna.txt", 
                                 organism=ORGANISM, binsize=BinSizes["bamCoverage"]),
        bamCoverage2_rna = expand("logs/metrics_bamCoverage2_{organism}_bs{binsize}_rna.txt", 
                                 organism=ORGANISM, binsize=BinSizes["bamCoverage"]),
        bamCoverage1_wgs = expand("logs/metrics_bamCoverage1_{organism}_bs{binsize}_wgs.txt", 
                                 organism=ORGANISM, binsize=BinSizes["bamCoverage"]),
        bamCoverage2_wgs = expand("logs/metrics_bamCoverage2_{organism}_bs{binsize}_wgs.txt", 
                                 organism=ORGANISM, binsize=BinSizes["bamCoverage"]),
        bamCompare1 = expand("logs/metrics_bamCompare1_{organism}_bs{binsize}.txt", 
                            organism=ORGANISM, binsize=BinSizes["bamCompare"]),
        bamCompare2 = expand("logs/metrics_bamCompare2_{organism}_bs{binsize}.txt", 
                            organism=ORGANISM, binsize=BinSizes["bamCompare"]),
        computeMatrix1 = expand("logs/metrics_computeMatrix1_{organism}_bs{binsize}.txt", 
                               organism=ORGANISM, binsize=BinSizes["computeMatrix"]),
        computeMatrix2 = expand("logs/metrics_computeMatrix2_{organism}_bs{binsize}.txt", 
                               organism=ORGANISM, binsize=BinSizes["computeMatrix"]),
        multiBamSummary1 = expand("logs/metrics_multiBamSummary1_{organism}_bs{binsize}.txt", 
                                 organism=ORGANISM, binsize=BinSizes["multiBamSummary"]),
        multiBamSummary2 = expand("logs/metrics_multiBamSummary2_{organism}_bs{binsize}.txt", 
                                 organism=ORGANISM, binsize=BinSizes["multiBamSummary"])
    output:
        # CPU Time plots
        bamCoverage_cputime_plot = f"output/bamCoverage_{ORGANISM}_bs{BinSizes['bamCoverage']}_cputime.png",
        bamCompare_cputime_plot = f"output/bamCompare_{ORGANISM}_bs{BinSizes['bamCompare']}_cputime.png",
        computeMatrix_cputime_plot = f"output/computeMatrix_{ORGANISM}_bs{BinSizes['computeMatrix']}_cputime.png",
        multiBamSummary_cputime_plot = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes['multiBamSummary']}_cputime.png",
        
        # Wall Time plots
        bamCoverage_wall_time_plot = f"output/bamCoverage_{ORGANISM}_bs{BinSizes['bamCoverage']}_walltime.png",
        bamCompare_wall_time_plot = f"output/bamCompare_{ORGANISM}_bs{BinSizes['bamCompare']}_walltime.png",
        computeMatrix_wall_time_plot = f"output/computeMatrix_{ORGANISM}_bs{BinSizes['computeMatrix']}_walltime.png",
        multiBamSummary_wall_time_plot = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes['multiBamSummary']}_walltime.png",
        
        # Memory plots
        bamCoverage_mem_plot = f"output/bamCoverage_{ORGANISM}_bs{BinSizes['bamCoverage']}_memory.png",
        bamCompare_mem_plot = f"output/bamCompare_{ORGANISM}_bs{BinSizes['bamCompare']}_memory.png",
        computeMatrix_mem_plot = f"output/computeMatrix_{ORGANISM}_bs{BinSizes['computeMatrix']}_memory.png",
        multiBamSummary_mem_plot = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes['multiBamSummary']}_memory.png"
    threads: 1
    params:
        # Pre-compute these paths to avoid string formatting issues in the shell command
        bamCoverage_output = f"output/bamCoverage_{ORGANISM}_bs{BinSizes['bamCoverage']}",
        bamCompare_output = f"output/bamCompare_{ORGANISM}_bs{BinSizes['bamCompare']}",
        computeMatrix_output = f"output/computeMatrix_{ORGANISM}_bs{BinSizes['computeMatrix']}",
        multiBamSummary_output = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes['multiBamSummary']}"
    shell:
        """
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.bamCoverage_output}.png \
            {input.bamCoverage1_chip},{input.bamCoverage1_rna},{input.bamCoverage1_wgs} \
            {input.bamCoverage2_chip},{input.bamCoverage2_rna},{input.bamCoverage2_wgs}
            
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.bamCompare_output}.png \
            {input.bamCompare1} {input.bamCompare2}
            
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.computeMatrix_output}.png \
            {input.computeMatrix1} {input.computeMatrix2}
            
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.multiBamSummary_output}.png \
            {input.multiBamSummary1} {input.multiBamSummary2}
        """

rule create_md:
    input:
        expand("output/bamCoverage_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["bamCoverage"], type=["walltime", "cputime", "memory"]),
        expand("output/bamCompare_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["bamCompare"], type=["walltime", "cputime", "memory"]),
        expand("output/computeMatrix_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["computeMatrix"], type=["walltime", "cputime", "memory"]),
        expand("output/multiBamSummary_{organism}_bs{binsize}_{type}.png",
               organism=ORGANISM, binsize=BinSizes["multiBamSummary"], type=["walltime", "cputime", "memory"])
    output:
        report = "output/report.md"
    params:
        organism = ORGANISM.capitalize(),
        threads = Nthreads,
        ntimes = Ntimes,
        bamCoverage_bs = BinSizes["bamCoverage"],
        bamCompare_bs = BinSizes["bamCompare"],
        computeMatrix_bs = BinSizes["computeMatrix"],
        multiBamSummary_bs = BinSizes["multiBamSummary"]
    shell:
        """
        cat > {output.report} << EOL

<style>
img {{
  max-width: 100%;
  height: auto;
}}
</style>

# {params.organism}

> Threads: {params.threads}
>
> Ntimes: {params.ntimes}
>

## bamCoverage

> Bin Size: {params.bamCoverage_bs}

### Walltime

![](./bamCoverage_{ORGANISM}_bs{params.bamCoverage_bs}_walltime.png)

### CPU Time

![](./bamCoverage_{ORGANISM}_bs{params.bamCoverage_bs}_cputime.png)

### Memory

![](./bamCoverage_{ORGANISM}_bs{params.bamCoverage_bs}_memory.png)

## bamCompare

> Bin Size: {params.bamCompare_bs}

### Walltime

![](./bamCompare_{ORGANISM}_bs{params.bamCompare_bs}_walltime.png)

### CPU Time

![](./bamCompare_{ORGANISM}_bs{params.bamCompare_bs}_cputime.png)

### Memory

![](./bamCompare_{ORGANISM}_bs{params.bamCompare_bs}_memory.png)

## computeMatrix

> Bin Size: {params.computeMatrix_bs}

### Walltime

![](./computeMatrix_{ORGANISM}_bs{params.computeMatrix_bs}_walltime.png)

### CPU Time

![](./computeMatrix_{ORGANISM}_bs{params.computeMatrix_bs}_cputime.png)

### Memory

![](./computeMatrix_{ORGANISM}_bs{params.computeMatrix_bs}_memory.png)

## multiBamSummary

> Bin Size: {params.multiBamSummary_bs}

### Walltime

![](./multiBamSummary_{ORGANISM}_bs{params.multiBamSummary_bs}_walltime.png)

### CPU Time

![](./multiBamSummary_{ORGANISM}_bs{params.multiBamSummary_bs}_cputime.png)

### Memory

![](./multiBamSummary_{ORGANISM}_bs{params.multiBamSummary_bs}_memory.png)
EOL
        """

rule render_md:
    input: "output/report.md"
    output: "output/report.html"
    conda: "extras.env.yaml"
    shell: "pandoc {input} -o {output}"
