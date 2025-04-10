#!/usr/bin/env python
# USAGE: #Generate_central_csv.py --folders /Volumes/ExtremeSSD/deeptools25/20250408_mac_bs5000_human:8,/Volumes/ExtremeSSD/deeptools25/mac_20250407_bs10_triticum:12 --output results
#python3 Generate_central_csv.py --folders /Volumes/ExtremeSSD/deeptools25/20250409_deep22_bs5000_triticum:64 --output results --system linux


import argparse
import glob
import os
import re
import warnings
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import logging
import platform

# Initialize logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Command patterns to look for
COMMAND_PATTERNS = {
    "bamCoverage": r"bamCoverage(\d+)_(\w+)_(\d+)\.txt",  # Updated pattern to match bamCoverage1_wgs_3.txt
    "bamCompare": r"bamCompare(\d+)_(\d+)\.txt",
    "computeMatrix": r"computeMatrix(\d+)_(\d+)\.txt",
    "multiBamSummary": r"multiBamSummary(\d+)_(\d+)\.txt",
}

def get_command_key(file_info):
    """Create a unique key for each command configuration"""
    if file_info["protocol"]:
        return f"{file_info['command']}{file_info['backend']}_{file_info['binsize']}_{file_info['protocol']}"
    else:
        return f"{file_info['command']}{file_info['backend']}_{file_info['binsize']}"

def expand_path(path):
    """Expand user home directory in path"""
    return os.path.expanduser(path)

def parse_file_path(filepath):
    """Extract command, backend, binsize, and protocol from file path"""
    basename = os.path.basename(filepath)

    for cmd, pattern in COMMAND_PATTERNS.items():
        match = re.match(pattern, basename)
        if match:
            backend = match.group(1)
            binsize = match.group(2)
            protocol = match.group(3) if len(match.groups()) > 2 else None
            return {
                "command": cmd,
                "backend": backend,
                "binsize": binsize,
                "protocol": protocol,
                "filename": basename,
            }

    return None

def parse_filename(filepath):
    parts = filepath.split("/")
    organism = parts[-2]  # Extract the second-to-last part of the path as the organism
    org=organism.split("_")[-1]
    filename = parts[-1]
    folder_name = parts[-2]
    mode = re.match(r"(.*)1", filename.split("_")[0]).group(1)  # Extract string before "1"

    if mode == "computeMatrix":
        mode = "computematrix"
        file_type = "gtf_4" 
    elif mode == "multiBamSummary":
        mode = "mbs"
        file_type = "bins_" + re.search(r"bs(\d+)", folder_name).group(1)
    elif mode == "bamCoverage":
        file_type = filename.split("_")[1] 
    else:
        file_type = "chip"

    return mode, org, file_type

def collect_files_by_command(folders):
    """
    Collect and organize files by command type across all CPU configs.

    Returns a dictionary mapping command types to lists of files from each CPU config.
    """
    # Dictionary to store files for each command across different CPU configurations
    # Structure: {command_key: {cpu_count: filepath}}
    command_files = {}

    for folder_info in folders:
        folder_path = folder_info["path"]
        cpu_count = folder_info["cpu_count"]

        # Find all files in the folder with "1" in the file title as the first number
        txt_files = glob.glob(os.path.join(folder_path, "*1_*.txt"))
        logger.info(
            f"Found {len(txt_files)} files with '1' in the title in {folder_path} (CPU count: {cpu_count})"
        )

        for txt_file in txt_files:
            file_info = parse_file_path(txt_file)
            if not file_info:
                logger.warning(f"Could not parse file name: {txt_file}")
                continue

            cmd_key = get_command_key(file_info)

            if cmd_key not in command_files:
                command_files[cmd_key] = {}

            command_files[cmd_key][cpu_count] = txt_file

    return command_files

def collect_time_data(folders, output_dir, system):
    """
    Collect time data from all log files in the given folders.
    """
    # First organize files by command across CPU configurations
    command_files = collect_files_by_command(folders)
    #open the log file to log the agregated data
    # Now process each command group    
    for folder_info in folders:
        folder_path = folder_info["path"]
        folder_name = os.path.basename(folder_path.rstrip("/"))
        file_path = os.path.join(output_dir, f"{folder_name}_benchmark_results.csv")
        with open(file_path, 'w') as file1:  # Changed from 'a' to 'w'
            file1.write("platform,mode,organism,type,threads,dt3,dt4\n")
            if system == "linux":
                system = "rhel8.8:x86_64"
            else:
                system = platform.system()
            print(f"System: {system}")
            for cmd_key, cpu_files in command_files.items():           
                print(f"Processing command: {cmd_key}")
                print(f"CPU files: {cpu_files}")
                for key, value in cpu_files.items():
                    #extract some info from the file name
                    mode, organism, file_type = parse_filename(value)
                    formatted_time_dt3=""
                    formatted_time_dt4=""
                    with open(value, 'r') as file:
                        for line in file:
                            if system == "Darwin":
                                if "real" in line:
                                    real_time_seconds = float(line.split()[0])
                                    minutes = int(real_time_seconds // 60)
                                    seconds = real_time_seconds % 60
                                    formatted_time_dt3 = f"{minutes}m{seconds:.3f}s"
                                    print(f"Formatted real time (dt3): {formatted_time_dt3}")
                                    break
                            else:
                                if "User time (seconds): " in line:
                                    real_time_seconds = float(line.split(": ")[1])
                                    minutes = int(real_time_seconds // 60)
                                    seconds = real_time_seconds % 60
                                    formatted_time_dt3 = f"{minutes}m{seconds:.3f}s"
                                    print(f"Formatted real time (dt3): {formatted_time_dt3}")
                                    break

                    second_file = value.replace("1_", "2_")
                    if os.path.exists(second_file):
                        with open(second_file, 'r') as file:
                            for line in file:
                                if system == "Darwin":
                                    if "real" in line:
                                        real_time_seconds = float(line.split()[0])
                                        minutes = int(real_time_seconds // 60)
                                        seconds = real_time_seconds % 60
                                        formatted_time_dt4 = f"{minutes}m{seconds:.3f}s"
                                        print(f"Formatted real time (dt4): {formatted_time_dt4}")
                                        break
                                else:
                                    if "User time (seconds): " in line:
                                        real_time_seconds = float(line.split(": ")[1])
                                        minutes = int(real_time_seconds // 60)
                                        seconds = real_time_seconds % 60
                                        formatted_time_dt4 = f"{minutes}m{seconds:.3f}s"
                                        print(f"Formatted real time (dt4): {formatted_time_dt4}")
                                        break
                    file1.write(f"{system.lower()},{mode.lower()},{organism.lower()},{file_type.lower()},{key},{formatted_time_dt3},{formatted_time_dt4}\n")

    return None

def main():

    parser = argparse.ArgumentParser(
        description="Compare memory usage across multiple CPU configurations"
    )
    parser.add_argument(
        "--folders",
        required=True,
        help="Comma-separated list of folders containing CSV files, each with its CPU count. "
        "Format: path1:cpucount1,path2:cpucount2,...",
    )
    parser.add_argument(
        "--output",
        default="memory_comparison",
        help="Output directory for plots and summary",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )

    parser.add_argument(
        "--system", "-s", type=str, help="Specify the system (e.g., linux, mac, windows)"
    )

    args = parser.parse_args()
    folders = []
    system = args.system if args.system else platform.system()
    for folder_spec in args.folders.split(","):
        parts = folder_spec.split(":")
        if len(parts) != 2:
            logger.error(
                f"Invalid folder specification: {folder_spec}. Expected format: path:cpucount"
            )
            return 1

        folder_path, cpu_count = parts
        # Expand user home directory in path
        folder_path = expand_path(folder_path)

        folders.append({"path": folder_path, "cpu_count": int(cpu_count)})

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Collect memory data
    logger.info("Collecting data from CSV files...")
    time_data = collect_time_data(folders, args.output, system)


if __name__ == "__main__":
    exit(main())

