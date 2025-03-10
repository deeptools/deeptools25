# Adjust these if you want
ORGANISM = "homo"
GTF = "regions/homo.v91.sample25k.gtf"

BinSizes = {
    "bamCoverage": 10,
    "bamCompare": 100,
    "computeMatrix": 300,
    "multiBamSummary": 5000,
}

Ntimes = 10
Nthreads = 16


# Do not edit any further
import platform
system = platform.system()
if system == "Linux":
    timeCmd = "/usr/bin/time -v"
elif system == "Darwin":
    timeCmd = "/usr/bin/time -al"
else:
    raise ValueError(f"Unsupported platform: {system}")

# Helper fn. to keep logs of failed jobs and verify file existence
shell.prefix("""
function on_error() {{ 
    cp $1 $1.failed.$(date +%Y%m%d_%H%M%S)
    echo "ERROR: Command failed, log saved as $1.failed.$(date +%Y%m%d_%H%M%S)" >&2
    return 1  # This will propagate the error
}}; 

function verify_file() {{
    local timeout=300  # Default timeout 5 minutes
    
    # Check if first arg is a number (timeout)
    if [[ $1 =~ ^[0-9]+$ ]]; then
        timeout=$1
        shift  # Remove first argument
    fi
    
    # Check if there are no files to check
    if [ $# -eq 0 ]; then
        echo "ERROR: No files specified for verification" >&2
        return 1
    fi
    
    local start_time=$(date +%s)
    local end_time=$((start_time + timeout))
    
    # Loop until timeout
    while [ $(date +%s) -lt $end_time ]; do
        local all_exist=true
        
        # Check each file
        for file in "$@"; do
            if [ ! -f "$file" ]; then
                all_exist=false
                echo "Waiting for file: $file" >&2
                break
            fi
        done
        
        # If all files exist, we're good
        if $all_exist; then
            return 0
        fi
        
        # Wait before checking again
        sleep 10
    done
    
    # Timeout occurred - report missing files
    echo "ERROR: Not all input files found after $timeout seconds:" >&2
    for file in "$@"; do
        if [ ! -f "$file" ]; then
            echo "  - Missing: $file" >&2
        fi
    done
    return 1
}};
"""
)

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
               organism=ORGANISM, binsize=BinSizes["bamCoverage"], type=["walltime", "cputime", "memory"]),
        expand("output/bamCompare_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["bamCompare"], type=["walltime", "cputime", "memory"]),
        expand("output/computeMatrix_{organism}_bs{binsize}_{type}.png", 
               organism=ORGANISM, binsize=BinSizes["computeMatrix"], type=["walltime", "cputime", "memory"]),
        expand("output/multiBamSummary_{organism}_bs{binsize}_{type}.png",
               organism=ORGANISM, binsize=BinSizes["multiBamSummary"], type=["walltime", "cputime", "memory"])


rule bamCoverage2:
    input:
        bam = lambda wildcards: FILES[ORGANISM + "_" + wildcards.protocol]
    output:
        bed = directory("output/bamCoverage2_{protocol}"),
        iter_file = "output/benchmark_iteration_bamCoverage2_{protocol}.txt",
        done = expand("output/bamCoverage2_{{protocol}}_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes['bamCoverage']}_{{protocol}}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCoverage"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        verify_file {input.bam} 600 || {{ echo "Timeout! Not all input files were found." 0 || exit 10 || exit 1 exit 1; }}
        mkdir -p {output.bed}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_file="{output.bed}/iter_${{curr_iter}}.bg"
        log_file="logs/bamCoverage2_{wildcards.protocol}_${{curr_iter}}.txt"
        
        {timeCmd} bamCoverage -b {input.bam} -o $out_file -of bedgraph \
            -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
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
        iter_file = "output/benchmark_iteration_bamCoverage1_{protocol}.txt",
        done = expand("output/bamCoverage1_{{protocol}}_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes['bamCoverage']}_{{protocol}}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCoverage"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        verify_file {input.bam} 600 || {{ echo "Timeout! Not all input files were found." 0 || exit 10 || exit 1 exit 1; }}
        mkdir -p {output.bed}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_file="{output.bed}/iter_${{curr_iter}}.bg"
        log_file="logs/bamCoverage1_{wildcards.protocol}_${{curr_iter}}.txt"
        
        {timeCmd} bamCoverage -b {input.bam} -o $out_file -of bedgraph \
            -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
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
        iter_file = "output/benchmark_iteration_bamCompare2.txt",
        done = expand("output/bamCompare2_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/bamCompare2_{ORGANISM}_bs{BinSizes['bamCompare']}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCompare"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        verify_file 600 {input.bam1} {input.bam2} || exit 1
        mkdir -p {output.bw}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_file="{output.bw}/iter_${{curr_iter}}.bw"
        log_file="logs/bamCompare2_${{curr_iter}}.txt"
        
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o $out_file -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
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
        iter_file = "output/benchmark_iteration_bamCompare1.txt",
        done = expand("output/bamCompare1_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/bamCompare1_{ORGANISM}_bs{BinSizes['bamCompare']}.txt", Ntimes)
    params:
        binsize = BinSizes["bamCompare"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        verify_file 1200 {input.bam1} {input.bam2} || exit 1
        mkdir -p {output.bw}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_file="{output.bw}/iter_${{curr_iter}}.bw"
        log_file="logs/bamCompare1_${{curr_iter}}.txt"
        
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o $out_file -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
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
        iter_file = "output/benchmark_iteration_computeMatrix2.txt",
        done = expand("output/computeMatrix2_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/computeMatrix2_{ORGANISM}_bs{BinSizes['computeMatrix']}.txt", Ntimes)
    params:
        binsize = BinSizes["computeMatrix"],
        upstream = 2 * BinSizes["computeMatrix"],
        downstream = 2 * BinSizes["computeMatrix"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        mkdir -p {output.npz}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_file="{output.npz}/iter_${{curr_iter}}.npz"
        log_file="logs/computeMatrix2_${{curr_iter}}.txt"
        
        # Try to match iterations, otherwise use first file
        if [ -f "{input.bw_dir}/iter_${{curr_iter}}.bw" ]; then
            input_bw="{input.bw_dir}/iter_${{curr_iter}}.bw"
        else
            input_bw=$(ls {input.bw_dir}/iter_*.bw | head -n 1)
        fi
        verify_file $input_bw 600 || {{ echo "Timeout! Not all input files were found." 0 || exit 10 || exit 1 exit 1; }}
        
        if [ "{ORGANISM}" = "homo" ]; then
          verify_file {humanBigwigs} 600 || {{ echo "Timeout! Not all input files were found." 0 || exit 10 || exit 1 exit 1; }}
          {timeCmd} computeMatrix reference-point \
              -S $input_bw {humanBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        else
          {timeCmd} computeMatrix reference-point \
              -S $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
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
        iter_file = "output/benchmark_iteration_computeMatrix1.txt",
        done = expand("output/computeMatrix1_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/computeMatrix1_{ORGANISM}_bs{BinSizes['computeMatrix']}.txt", Ntimes)
    params:
        binsize = BinSizes["computeMatrix"],
        upstream = 2 * BinSizes["computeMatrix"],
        downstream = 2 * BinSizes["computeMatrix"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        mkdir -p {output.npz}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_file="{output.npz}/iter_${{curr_iter}}.npz"
        log_file="logs/computeMatrix1_${{curr_iter}}.txt"
        
        # Try to match iterations, otherwise use first file
        if [ -f "{input.bw_dir}/iter_${{curr_iter}}.bw" ]; then
            input_bw="{input.bw_dir}/iter_${{curr_iter}}.bw"
        else
            input_bw=$(ls {input.bw_dir}/iter_*.bw | head -n 1)
        fi
        verify_file $input_bw 600 || {{ echo "Timeout! Not all input files were found." 0 || exit 10 || exit 1 exit 1; }}
        
        if [ "{ORGANISM}" = "homo" ]; then
          verify_file {humanBigwigs} 600 || {{ echo "Timeout! Not all input files were found." 0 || exit 10 || exit 1 exit 1; }}
          {timeCmd} computeMatrix reference-point \
              -S $input_bw {humanBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        else
          {timeCmd} computeMatrix reference-point \
              -S $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
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
        iter_file = "output/benchmark_iteration_multiBamSummary2.txt",
        done = expand("output/multiBamSummary2_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/multiBamSummary2_{ORGANISM}_bs{BinSizes['multiBamSummary']}.txt", Ntimes)
    params:
        binsize = BinSizes["multiBamSummary"]
    threads: Nthreads
    conda: "v4.env.yaml"
    shell:
        """
        verify_file 600 {input.bam} || exit 1
        mkdir -p {output.npz}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_npz="{output.npz}/iter_${{curr_iter}}.npz"
        out_raw="{output.npz}/iter_${{curr_iter}}.outraw.tab"
        log_file="logs/multiBamSummary2_${{curr_iter}}.txt"
        
        if [ "{ORGANISM}" = "homo" ]; then
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        else
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
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
        iter_file = "output/benchmark_iteration_multiBamSummary1.txt",
        done = expand("output/multiBamSummary1_done_{iter}.txt", iter=range(1, Ntimes+1))
    benchmark:
        repeat(f"logs/multiBamSummary1_{ORGANISM}_bs{BinSizes['multiBamSummary']}.txt", Ntimes)
    params:
        binsize = BinSizes["multiBamSummary"]
    threads: Nthreads
    conda: "v3.env.yaml"
    shell:
        """
        verify_file 600 {input.bam} || exit 1
        mkdir -p {output.npz}
        curr_iter=$(cat {output.iter_file} 2>/dev/null || echo 1)
        out_npz="{output.npz}/iter_${{curr_iter}}.npz"
        out_raw="{output.npz}/iter_${{curr_iter}}.outraw.tab"
        log_file="logs/multiBamSummary1_${{curr_iter}}.txt"
        
        if [ "{ORGANISM}" = "homo" ]; then
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        else
            {timeCmd} multiBamSummary bins -b {input.bam} {input.bam} {input.bam} \
                -o $out_npz --outRawCounts $out_raw -bs {params.binsize} -p {threads} \
                    > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
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
        
        # Original input files
        bamCoverage1_chip = f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes['bamCoverage']}_chip.txt",
        bamCoverage2_chip = f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes['bamCoverage']}_chip.txt",
        bamCoverage1_rna = f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes['bamCoverage']}_rna.txt",
        bamCoverage2_rna = f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes['bamCoverage']}_rna.txt",
        bamCoverage1_wgs = f"logs/bamCoverage1_{ORGANISM}_bs{BinSizes['bamCoverage']}_wgs.txt",
        bamCoverage2_wgs = f"logs/bamCoverage2_{ORGANISM}_bs{BinSizes['bamCoverage']}_wgs.txt",
        bamCompare1 = f"logs/bamCompare1_{ORGANISM}_bs{BinSizes['bamCompare']}.txt",
        bamCompare2 = f"logs/bamCompare2_{ORGANISM}_bs{BinSizes['bamCompare']}.txt",
        computeMatrix1 = f"logs/computeMatrix1_{ORGANISM}_bs{BinSizes['computeMatrix']}.txt",
        computeMatrix2 = f"logs/computeMatrix2_{ORGANISM}_bs{BinSizes['computeMatrix']}.txt",
        multiBamSummary1 = f"logs/multiBamSummary1_{ORGANISM}_bs{BinSizes['multiBamSummary']}.txt",
        multiBamSummary2 = f"logs/multiBamSummary2_{ORGANISM}_bs{BinSizes['multiBamSummary']}.txt"
    output:
        # CPU Time plots
        bamCoverage_cpu_time_plot = f"output/bamCoverage_{ORGANISM}_bs{BinSizes['bamCoverage']}_cputime.png",
        bamCompare_cpu_time_plot = f"output/bamCompare_{ORGANISM}_bs{BinSizes['bamCompare']}_cputime.png",
        computeMatrix_cpu_time_plot = f"output/computeMatrix_{ORGANISM}_bs{BinSizes['computeMatrix']}_cputime.png",
        multiBamSummary_cpu_time_plot = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes['multiBamSummary']}_cputime.png",
        
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
    params:
        # Pre-compute these paths to avoid string formatting issues in the shell command
        bamCoverage_output = f"output/bamCoverage_{ORGANISM}_bs{BinSizes['bamCoverage']}",
        bamCompare_output = f"output/bamCompare_{ORGANISM}_bs{BinSizes['bamCompare']}",
        computeMatrix_output = f"output/computeMatrix_{ORGANISM}_bs{BinSizes['computeMatrix']}",
        multiBamSummary_output = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes['multiBamSummary']}"
    shell:
        """
        verify_file {input} 600 || {{ echo "Timeout! Not all input files were found." 0 || exit 10 || exit 1 exit 1; }}
        python3 present_results.py --ntimes {Ntimes} {params.bamCoverage_output}.png \
            {input.bamCoverage1_chip},{input.bamCoverage1_rna},{input.bamCoverage1_wgs} \
            {input.bamCoverage2_chip},{input.bamCoverage2_rna},{input.bamCoverage2_wgs}
        python3 present_results.py --ntimes {Ntimes} {params.bamCompare_output}.png {input.bamCompare1} {input.bamCompare2}
        python3 present_results.py --ntimes {Ntimes} {params.computeMatrix_output}.png {input.computeMatrix1} {input.computeMatrix2}
        python3 present_results.py --ntimes {Ntimes} {params.multiBamSummary_output}.png {input.multiBamSummary1} {input.multiBamSummary2}
        """
