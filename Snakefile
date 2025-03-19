BinSizes = {
    "bamCoverage": 10,
    "bamCompare": 100,
    "computeMatrix": 300,
    "multiBamSummary": 5000,
}

Ntimes = 10
Nthreads = 64

ORGANISM = config.get("organism", "homo")

# Do not edit any further
if ORGANISM == "homo":
    GTF = f"regions/{ORGANISM}.v91.sample25k.gtf"
elif ORGANISM == "triticum":
    GTF = f"regions/{ORGANISM}.v60.sample25k.gtf"
else:
    raise ValueError(f"Unsupported organism: {ORGANISM}")

# Apply resource adjustments based on organism, but only when running via slurm
def set_organism_resources():
    """Apply organism-specific resource configurations from config.yaml when using Slurm executor"""
    import yaml
    import os
    
    # Check if we're likely running with the Slurm executor
    is_slurm = False
    
    # Method 1: Check if we're running under a Slurm profile
    if workflow.config.get("cluster", None) or workflow.config.get("profile", "").endswith("slurm"):
        is_slurm = True
    
    # Method 2: Check for common Snakemake-to-Slurm environment variables
    if os.environ.get("SNAKEMAKE_CLUSTER_SLURM") or os.environ.get("SNAKEMAKE_PROFILE", "").endswith("slurm"):
        is_slurm = True
        
    # Method 3: Check for specific Slurm variables that might be present
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CLUSTER_NAME"):
        is_slurm = True
    
    if not is_slurm:
        print("Not running via Slurm executor - skipping organism-specific resource configuration")
        return
        
    print("Running via Slurm executor - applying organism-specific resources")
    
    # First, get the organism-specific resources
    organism_resources = None
    try:
        config_path = "snk-slurm-exe/config.yaml"
        if not os.path.exists(config_path):
            print(f"Warning: Config file {config_path} not found, skipping organism resources")
            return
            
        with open(config_path, "r") as f:
            config_data = yaml.safe_load(f)
            if "organism_resources" in config_data and ORGANISM in config_data["organism_resources"]:
                organism_resources = config_data["organism_resources"][ORGANISM]
                print(f"Found {ORGANISM}-specific resources in config")
            else:
                print(f"No organism-specific resources found for {ORGANISM}")
                return
    except Exception as e:
        print(f"Warning: Could not load organism-specific resources: {e}")
        return
    
    # If organism-specific resources exist, update rule resources
    if organism_resources:
        for rule_name, rule_resources in organism_resources.items():
            if rule_name in workflow.rules:
                snakemake_rule = workflow.rules[rule_name]
                for resource_name, value in rule_resources.items():
                    snakemake_rule.resources[resource_name] = value
                print(f"Updated resources for {rule_name} with {ORGANISM}-specific settings: {rule_resources}")

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

function verify_file() {{
    local timeout=300  # Default timeout 5 minutes
    local verbose=false
    
    # Process options
    while [[ "$1" == "-"* ]]; do
        case "$1" in
            -t|--timeout)
                timeout=$2
                shift 2
                ;;
            -v|--verbose)
                verbose=true
                shift
                ;;
            *)
                echo "Unknown option: $1" >&2
                return 1
                ;;
        esac
    done
    
    # Check if there are no files to check
    if [ $# -eq 0 ]; then
        echo "ERROR: No files specified for verification" >&2
        return 1
    fi
    
    $verbose && echo "Verifying $(($#)) files with timeout of $timeout seconds..."
    
    local start_time=$(date +%s)
    local end_time=$((start_time + timeout))
    local check_count=0
    
    # First quick check - if all files exist right away, return immediately
    local missing=0
    for file in "$@"; do
        if [ ! -f "$file" ] || [ ! -r "$file" ]; then
            missing=$((missing + 1))
            $verbose && echo "Initially missing: $file"
        fi
    done
    
    if [ $missing -eq 0 ]; then
        $verbose && echo "All files present immediately."
        return 0
    fi
    
    # Loop until timeout
    while [ $(date +%s) -lt $end_time ]; do
        check_count=$((check_count + 1))
        missing=0
        
        # Check each file
        for file in "$@"; do
            if [ ! -f "$file" ] || [ ! -r "$file" ]; then
                missing=$((missing + 1))
                if $verbose && [ $((check_count % 6)) -eq 0 ]; then
                    # Only print every ~60 seconds if verbose
                    echo "Still waiting for: $file"
                fi
            fi
        done
        
        # If all files exist and are readable, we're good
        if [ $missing -eq 0 ]; then
            elapsed=$(($(date +%s) - start_time))
            echo "All files present after $elapsed seconds."
            return 0
        fi
        
        # Wait before checking again (but only print status occasionally)
        if $verbose && [ $((check_count % 6)) -eq 0 ]; then
            elapsed=$(($(date +%s) - start_time))
            remaining=$((timeout - elapsed))
            echo "Still missing $missing files after $elapsed seconds. Will wait $remaining more seconds."
        fi
        
        sleep 10
    done
    
    # Timeout occurred - report missing files
    echo "ERROR: Not all input files found after $timeout seconds:" >&2
    for file in "$@"; do
        if [ ! -f "$file" ] || [ ! -r "$file" ]; then
            echo "  - Missing or unreadable: $file" >&2
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
        verify_file -t 600 {input.bam} || {{ echo "Timeout! Not all input files were found." && exit 1; }}
        mkdir -p {output.bed}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bed}/iter_${{curr_iter}}.bg"
        log_file="logs/bamCoverage2_{wildcards.protocol}_${{curr_iter}}.txt"
        
        {timeCmd} bamCoverage -b {input.bam} -o $out_file -of bedgraph \
            -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
        # Verify output files were created properly
        verify_file -t 3600 $out_file || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
        
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
        verify_file -t 600 {input.bam} || {{ echo "Timeout! Not all input files were found." && exit 1; }}
        mkdir -p {output.bed}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bed}/iter_${{curr_iter}}.bg"
        log_file="logs/bamCoverage1_{wildcards.protocol}_${{curr_iter}}.txt"
        
        {timeCmd} bamCoverage -b {input.bam} -o $out_file -of bedgraph \
            -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
        # Verify output files were created properly
        verify_file -t 3600 $out_file || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
        
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
        verify_file -t 600 {input.bam1} {input.bam2} || exit 1
        mkdir -p {output.bw}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bw}/iter_${{curr_iter}}.bw"
        log_file="logs/bamCompare2_${{curr_iter}}.txt"
        
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o $out_file -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
        # Verify output files were created properly
        verify_file -t 3600 $out_file || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
                
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
        verify_file -t 1200 {input.bam1} {input.bam2} || exit 1
        mkdir -p {output.bw}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
        out_file="{output.bw}/iter_${{curr_iter}}.bw"
        log_file="logs/bamCompare1_${{curr_iter}}.txt"
        
        {timeCmd} bamCompare -b1 {input.bam1} -b2 {input.bam2} \
            -o $out_file -bs {params.binsize} -p {threads} \
                > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
                
        # Verify output files were created properly
        verify_file -t 3600 $out_file || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
                
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
          verify_file -t 600 $input_bw {humanBigwigs} || {{ echo "Timeout! Not all input files were found." && exit 1; }}
          {timeCmd} computeMatrix reference-point \
              -S $input_bw {humanBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        else
          verify_file -t 600 $input_bw || {{ echo "Timeout! Not all input files were found." && exit 1; }}
          {timeCmd} computeMatrix reference-point \
              -S $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        fi
        
        # Verify output files were created properly
        verify_file -t 3600 $out_file || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
        
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
          verify_file -t 600 $input_bw {humanBigwigs} || {{ echo "Timeout! Not all input files were found." && exit 1; }}
          {timeCmd} computeMatrix reference-point \
              -S $input_bw {humanBigwigs} \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        else
          verify_file -t 600 $input_bw || {{ echo "Timeout! Not all input files were found." && exit 1; }}
          {timeCmd} computeMatrix reference-point \
              -S $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw $input_bw \
              -R {input.gtf} -o $out_file -a {params.downstream} -b {params.upstream} -bs {params.binsize} -p {threads} --missingDataAsZero \
                  > $log_file 2>&1 || {{ on_error $log_file && exit 1; }}
        fi
        
        # Verify output files were created properly
        verify_file -t 3600 $out_file || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
        
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
        verify_file -t 600 {input.bam} || exit 1
        mkdir -p {output.npz}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
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
        
        # Verify output files were created properly
        verify_file -t 3600 $out_raw $out_npz || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
        
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
        verify_file -t 600 {input.bam} || exit 1
        mkdir -p {output.npz}
        [ ! -f {output.iter_file} ] && echo "1" > {output.iter_file} || :
        curr_iter=$(cat {output.iter_file})
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
        
        # Verify output files were created properly
        verify_file -t 3600 $out_raw $out_npz || {{ on_error $log_file && echo "Output file verification failed!" && exit 1; }}
        
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
    threads: 1
    params:
        # Pre-compute these paths to avoid string formatting issues in the shell command
        bamCoverage_output = f"output/bamCoverage_{ORGANISM}_bs{BinSizes['bamCoverage']}",
        bamCompare_output = f"output/bamCompare_{ORGANISM}_bs{BinSizes['bamCompare']}",
        computeMatrix_output = f"output/computeMatrix_{ORGANISM}_bs{BinSizes['computeMatrix']}",
        multiBamSummary_output = f"output/multiBamSummary_{ORGANISM}_bs{BinSizes['multiBamSummary']}"
    shell:
        """
        verify_file -t 600 {input} || {{ echo "Timeout! Not all input files were found." && exit 1; }}
        
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.bamCoverage_output}.png \
            {input.bamCoverage1_chip},{input.bamCoverage1_rna},{input.bamCoverage1_wgs} \
            {input.bamCoverage2_chip},{input.bamCoverage2_rna},{input.bamCoverage2_wgs}
            
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.bamCompare_output}.png \
            {input.bamCompare1} {input.bamCompare2}
            
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.computeMatrix_output}.png \
            {input.computeMatrix1} {input.computeMatrix2}
            
        python3 present_results.py --threads {Nthreads} --ntimes {Ntimes} {params.multiBamSummary_output}.png \
            {input.multiBamSummary1} {input.multiBamSummary2}
        
        # Verify all output files were created
        verify_file -t 300 \
            {output.bamCoverage_cpu_time_plot} {output.bamCoverage_wall_time_plot} {output.bamCoverage_mem_plot} \
            {output.bamCompare_cpu_time_plot} {output.bamCompare_wall_time_plot} {output.bamCompare_mem_plot} \
            {output.computeMatrix_cpu_time_plot} {output.computeMatrix_wall_time_plot} {output.computeMatrix_mem_plot} \
            {output.multiBamSummary_cpu_time_plot} {output.multiBamSummary_wall_time_plot} {output.multiBamSummary_mem_plot} \
            || {{ echo "Output plot verification failed!" && exit 1; }}
        """
