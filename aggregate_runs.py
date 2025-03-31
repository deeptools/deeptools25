#!/usr/bin/env python
# USAGE: $ aggregate_runs.py --folders ~/data/bioinfo/benchmarking/hpc12_20250325235126_homo:12,~/data/bioinfo/benchmarking/hpc16_20250325001553_homo:16 --output memory_comparison
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

# Matplotlib settings for our use-case(s)
mpl.use("Agg")

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from matplotlib.cbook import MatplotlibDeprecationWarning

    warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)
except ImportError:
    warnings.filterwarnings("ignore", message=".*deprecated.*", category=Warning)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger("memory-compare")

# Command patterns to look for
COMMAND_PATTERNS = {
    "bamCoverage": r"bamCoverage(\d+)_(\d+)_(\w+)\.csv",
    "bamCompare": r"bamCompare(\d+)_(\d+)\.csv",
    "computeMatrix": r"computeMatrix(\d+)_(\d+)\.csv",
    "multiBamSummary": r"multiBamSummary(\d+)_(\d+)\.csv",
}

# Colors for different CPU configurations
COLOR_CYCLE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#aec7e8",
    "#ffbb78",
]


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


def get_command_key(file_info):
    """Create a unique key for each command configuration"""
    if file_info["protocol"]:
        return f"{file_info['command']}{file_info['backend']}_{file_info['binsize']}_{file_info['protocol']}"
    else:
        return f"{file_info['command']}{file_info['backend']}_{file_info['binsize']}"


def read_memory_from_csv(filepath):
    """Read memory data from CSV file"""
    try:
        df = pd.read_csv(filepath)
        # Look for the memory column - might be kernel_max_rss or python_memory_mb
        memory_col = None
        if "kernel_max_rss" in df.columns:
            memory_col = "kernel_max_rss"
        elif "python_memory_mb" in df.columns:
            memory_col = "python_memory_mb"

        if memory_col is None:
            logger.warning(f"No memory column found in {filepath}")
            return []

        # Convert to list and filter out NaN/None values
        memory_values = df[memory_col].dropna().tolist()

        # Log the number of measurements
        logger.info(
            f"Found {len(memory_values)} memory measurements in {os.path.basename(filepath)}"
        )

        return memory_values
    except Exception as e:
        logger.error(f"Error reading {filepath}: {e}")
        return []


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


def collect_memory_data(folders):
    """
    Collect memory data from all CSV files in the given folders.
    Group by command type for comparison across CPU counts.
    """
    # First organize files by command across CPU configurations
    command_files = collect_files_by_command(folders)

    # Now process each command group
    memory_data = {}

    for cmd_key, cpu_files in command_files.items():
        # Only include commands that have data from at least 2 CPU configs
        if len(cpu_files) >= 2:
            memory_data[cmd_key] = {}

            for cpu_count, filepath in cpu_files.items():
                memory_values = read_memory_from_csv(filepath)
                if memory_values:
                    memory_data[cmd_key][cpu_count] = memory_values
                else:
                    logger.warning(f"No valid memory data in: {filepath}")
        else:
            logger.warning(
                f"Skipping {cmd_key}: data available for only {len(cpu_files)} CPU configuration(s)"
            )

    return memory_data


def plot_memory_comparison(cmd_key, memory_data, output_dir):
    """Create a boxplot comparing memory usage across different CPU counts"""
    cpu_counts = sorted(memory_data.keys())

    # Prepare plot data
    plot_data = [memory_data[cpu] for cpu in cpu_counts]
    labels = [f"{cpu} CPUs (n={len(memory_data[cpu])})" for cpu in cpu_counts]

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 7))

    bp = ax.boxplot(plot_data, patch_artist=True, labels=labels)

    # Add styling to boxplot
    for i, patch in enumerate(bp["boxes"]):
        color = COLOR_CYCLE[i % len(COLOR_CYCLE)]
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    plt.setp(bp["medians"], color="black")
    plt.setp(bp["whiskers"], color="navy")
    plt.setp(bp["caps"], color="navy")
    plt.setp(bp["fliers"], marker="o", markerfacecolor="red", alpha=0.5)

    # Add individual data points with jitter for better visibility
    for i, data in enumerate(plot_data):
        x = np.random.normal(i + 1, 0.04, len(data))
        ax.scatter(x, data, alpha=0.4, s=20, color="darkblue")

    # Add mean and median labels
    for i, data in enumerate(plot_data):
        mean_val = np.mean(data)
        median_val = np.median(data)
        y_pos = np.max(data) * 1.05
        ax.text(i + 1, y_pos, f"Mean: {mean_val:.1f}", ha="center")
        ax.text(i + 1, y_pos * 0.95, f"Median: {median_val:.1f}", ha="center")

    # Clean up command key for display
    display_name = cmd_key.replace("_", " ")
    match = re.match(r"(\w+)(\d+)_(.*)", cmd_key)
    if match:
        cmd, backend, rest = match.groups()
        if "_" in rest:
            parts = rest.split("_")
            binsize = parts[0]
            protocol = parts[1]
            display_name = (
                f"{cmd} (Backend {backend}, Bin Size {binsize}, {protocol.upper()})"
            )
        else:
            display_name = f"{cmd} (Backend {backend}, Bin Size {rest})"

    ax.set_title(f"Memory Usage Comparison for {display_name}")
    ax.set_ylabel("Memory (MB)")
    ax.set_xlabel("CPU Configuration")
    ax.grid(True, linestyle="--", alpha=0.7)

    # Add memory efficiency analysis
    memory_efficiency = analyze_memory_scaling(cpu_counts, plot_data)
    plt.figtext(
        0.5,
        0.01,
        memory_efficiency,
        ha="center",
        fontsize=10,
        bbox={"facecolor": "lightyellow", "alpha": 0.5, "pad": 5},
    )

    # Save plot
    output_file = os.path.join(output_dir, f"memory_comparison_{cmd_key}.png")
    plt.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)

    logger.info(f"Saved plot to {output_file}")

    return output_file


def analyze_memory_scaling(cpu_counts, data_sets):
    """Analyze memory scaling efficiency across CPU counts"""
    means = [np.mean(data) for data in data_sets]

    if len(means) < 2:
        return "Insufficient data for memory scaling analysis"

    # Calculate scaling factors
    base_cpu = cpu_counts[0]
    base_mem = means[0]

    scaling_text = []
    scaling_text.append("Memory Scaling Analysis:")

    for i in range(1, len(cpu_counts)):
        cpu = cpu_counts[i]
        mem = means[i]
        cpu_ratio = cpu / base_cpu
        mem_ratio = mem / base_mem
        efficiency = mem_ratio / cpu_ratio

        scaling_text.append(
            f"{base_cpu}→{cpu} CPUs: Memory {mem_ratio:.2f}x (CPU ratio: {cpu_ratio:.2f}x, "
            f"Efficiency: {efficiency:.2f})"
        )

        if efficiency < 0.8:
            scaling_text.append("👍 Good memory scaling (sub-linear)")
        elif efficiency < 1.1:
            scaling_text.append("🔄 Linear memory scaling")
        else:
            scaling_text.append("⚠️ Poor memory scaling (super-linear)")

    return "\n".join(scaling_text)


def generate_summary_report(results, output_file, used_folders):
    """Generate a summary report of all plots created"""
    with open(output_file, "w") as f:
        f.write("# Memory Usage Comparison Summary\n\n")
        f.write(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Add information about the folders used
        f.write("## CPU Configurations Analyzed\n\n")
        for folder in used_folders:
            f.write(f"- {folder['cpu_count']} CPUs: `{folder['path']}`\n")
        f.write("\n")

        f.write("## Tools Analyzed\n\n")

        # Group by command, protocol, and binsize
        by_command_protocol = {}
        for cmd_key in results:
            match = re.match(r"(\w+)(\d+)_(.*)", cmd_key)
            if match:
                cmd, backend, rest = match.groups()

                # Split rest into binsize and protocol if protocol exists
                if "_" in rest:
                    binsize, protocol = rest.split("_")
                    # Include protocol in the key to treat each protocol as separate
                    command_key = f"{cmd}_{protocol}"
                else:
                    binsize = rest
                    command_key = cmd

                # Add binsize to the key so different binsizes are grouped separately
                command_key = f"{command_key}_bs{binsize}"

                if command_key not in by_command_protocol:
                    by_command_protocol[command_key] = []
                by_command_protocol[command_key].append((cmd_key, results[cmd_key]))

        # Write details for each command
        for command_key, items in sorted(by_command_protocol.items()):
            # Format the section header more nicely
            if "_bs" in command_key:
                parts = command_key.split("_bs")
                cmd_part = parts[0].replace("_", " ")
                binsize = parts[1]
                header = f"{cmd_part} (Bin Size {binsize})"
            else:
                header = command_key.replace("_", " ")

            f.write(f"### {header}\n\n")

            for cmd_key, plot_path in sorted(items):
                match = re.match(r"(\w+)(\d+)_(.*)", cmd_key)
                if match:
                    cmd, backend, rest = match.groups()
                    if "_" in rest:
                        display_name = f"{cmd} (Backend {backend})"
                    else:
                        display_name = f"{cmd} (Backend {backend})"

                f.write(f"#### {display_name}\n\n")
                f.write(
                    f"![Memory comparison for {display_name}]({os.path.basename(plot_path)})\n\n"
                )

    logger.info(f"Generated summary report: {output_file}")


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

    # Set logging level based on verbose flag
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Parse folder paths and CPU counts
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
    logger.info("Collecting memory data from CSV files...")
    memory_data = collect_memory_data(folders)

    if not memory_data:
        logger.error("No valid comparison data found between the CPU configurations")
        return 1

    logger.info(f"Found data for {len(memory_data)} command configurations")

    # Generate plots for each command configuration
    results = {}
    for cmd_key, cpu_data in memory_data.items():
        if len(cpu_data) < 2:
            logger.warning(
                f"Skipping {cmd_key}: data available for fewer than 2 CPU configurations"
            )
            continue

        logger.info(f"Generating plot for {cmd_key}")
        output_file = plot_memory_comparison(cmd_key, cpu_data, args.output)
        results[cmd_key] = output_file

    # Generate summary report
    if results:
        summary_file = os.path.join(args.output, "memory_comparison_summary.md")
        generate_summary_report(results, summary_file, folders)

        logger.info(f"Completed analysis. Generated {len(results)} plots.")
        logger.info(f"Summary report: {summary_file}")
    else:
        logger.error(
            "No plots were generated. Check if there are matching files across CPU configurations."
        )
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
