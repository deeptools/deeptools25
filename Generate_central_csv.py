#!/usr/bin/env python
# USAGE: $ Generate_central_csv.py --folders ~/data/bioinfo/benchmarking/hpc12_20250325235126_homo:12,~/data/bioinfo/benchmarking/hpc16_20250325001553_homo:16 --output results
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


# Command patterns to look for
COMMAND_PATTERNS = {
    "bamCoverage": r"bamCoverage(\d+)_(\d+)_(\w+)\.csv",
    "bamCompare": r"bamCompare(\d+)_(\d+)\.csv",
    "computeMatrix": r"computeMatrix(\d+)_(\d+)\.csv",
    "multiBamSummary": r"multiBamSummary(\d+)_(\d+)\.csv",
}

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

        # Find all CSV files in the folder
        csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
        logger.info(
            f"Found {len(csv_files)} CSV files in {folder_path} (CPU count: {cpu_count})"
        )

        for csv_file in csv_files:
            file_info = parse_file_path(csv_file)
            if not file_info:
                logger.warning(f"Could not parse file name: {csv_file}")
                continue

            cmd_key = get_command_key(file_info)

            if cmd_key not in command_files:
                command_files[cmd_key] = {}

            command_files[cmd_key][cpu_count] = csv_file

    return command_files

def collect_time_data(folders):
    """
    Collect time data from all CSV files in the given folders.
    Group by command type for comparison across CPU counts.
    """
    # First organize files by command across CPU configurations
    command_files = collect_files_by_command(folders)
    #open the log file to log the agregated data
    # Now process each command group
    time_data = {}
    file_path= os.path.join(args.output, "BenchmarkResults.csv")
    with open(file_path, 'a') as file:
        file.write("platform,mode,organism,type,threads,dt3,dt4\n")
        system = platform.system()
        for cmd_key, cpu_files in command_files.items():
            # Only include commands that have data from at least 2 CPU configs
            if len(cpu_files) >= 2:
                memory_data[cmd_key] = {}

                for cpu_count, filepath in cpu_files.items():
                    memory_values = read_memory_from_csv(filepath)
                    if memory_values:
                        memory_data[cmd_key][cpu_count] = memory_values
                        file.write(
                            f"{system},{cmd_key},{cpu_files[cpu_count]},{cpu_count},{memory_values['dt3']},{memory_values['dt4']}\n"
                        )
                    else:
                        logger.warning(f"No valid memory data in: {filepath}")
            else:
                logger.warning(
                    f"Skipping {cmd_key}: data available for only {len(cpu_files)} CPU configuration(s)"
                )

    return memory_data

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

    args = parser.parse_args()
    folders = []
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
    time_data = collect_time_data(folders)


if __name__ == "__main__":
    exit(main())

