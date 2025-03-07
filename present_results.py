#!/usr/bin/env python
import argparse
import glob
import matplotlib.pyplot as plt
import pandas as pd
import platform
import re
import os
from pathlib import Path


def extract_metadata_from_path(file_path):
    try:
        match = re.search(
            r"/(\w+?)(\d+)_(?:.*?)(?:_bs(\d+))?(?:_(\w+))?\.txt", file_path
        )
        command = match.group(1)  # 'bamCompare' or 'bamCoverage'
        backend = match.group(2)  # '1'
        binsize = match.group(3)  # '100'
        protocol = match.group(4)  # 'chip' or None
        return command, backend, binsize, protocol
    except (AttributeError, IndexError) as e:
        raise ValueError(f"Could not parse filename pattern from {file_path}: {e}")


def parse_benchmark_file(file_path):
    times, memory, cpu_usage, io_in, io_out, mean_load, cpu_time, max_uss, max_pss = (
        [] for _ in range(9)
    )

    error_detected = False
    error_message = None

    # Read file just once for both error detection and parsing
    try:
        with open(file_path, "r") as file:
            content = file.readlines()

            # Check for global errors in the entire content
            full_content = "".join(content)
            if "does not exist" in full_content:
                error_detected = True
                error_message = "File not found error"
                print(f"Warning: {error_message} found in {file_path}")
            elif "Read block operation failed" in full_content:
                error_detected = True
                error_message = "Read block operation failed"
                print(f"Warning: {error_message} found in {file_path}")

            # Skip header line
            content = content[1:]

            # Parse lines (excluding error message lines)
            for line_num, line in enumerate(content, 2):
                # Skip error message lines
                if "does not exist" in line or "Read block operation failed" in line:
                    continue

                try:
                    parts = line.strip().split("\t")

                    # Skip lines with incomplete data
                    if len(parts) < 16:
                        continue

                    # For lines with zeros, add None values instead
                    if float(parts[0]) == 0 and float(parts[2]) == 0:
                        times.append(None)
                        memory.append(None)
                        max_uss.append(None)
                        max_pss.append(None)
                        io_in.append(None)
                        io_out.append(None)
                        mean_load.append(None)
                        cpu_time.append(None)
                        cpu_usage.append(None)
                        continue

                    times.append(float(parts[0]))
                    memory.append(float(parts[2]))
                    max_uss.append(float(parts[4]))
                    max_pss.append(float(parts[5]))
                    io_in.append(float(parts[6]))
                    io_out.append(float(parts[7]))
                    mean_load.append(float(parts[8]))
                    cpu_time.append(float(parts[9]))
                    cpu_usage.append(float(parts[15].split()[0]))
                except (IndexError, ValueError, TypeError) as e:
                    print(f"Warning: Error in line {line_num}: {e}")
                    continue
    except FileNotFoundError:
        error_detected = True
        error_message = "File not found"
        print(f"Warning: {error_message}: {file_path}")

    results = {
        "times": times,
        "memory": memory,
        "cpu_usage": cpu_usage,
        "io_in": io_in,
        "io_out": io_out,
        "mean_load": mean_load,
        "cpu_time": cpu_time,
        "max_uss": max_uss,
        "max_pss": max_pss,
    }

    # If we detected errors and got no data, add a None value
    if error_detected and not times:
        for key in results:
            results[key] = [None]
        results["error"] = [error_message]

    return results, error_detected


def handle_macos_memory(command, backend, protocol, memory):
    if not platform.system() == "Darwin":
        return memory

    Ntimes = len(memory)
    assert sum(memory) == 0, "Memory values are not all zeroes"

    if command == "bamCoverage":
        assert protocol is not None, (
            f"We couldn't retrieve protocol for {command}{backend}"
        )
        log_files = glob.glob(f"logs/{command}{backend}_{protocol}_[0-9]*.txt")
    else:
        log_files = glob.glob(f"logs/{command}{backend}_[0-9]*.txt")

    memory = parse_memory_from_logs(log_files)
    assert len(memory) == Ntimes, (
        f"Expected {Ntimes} memory values but got {len(memory)}, {memory}"
    )

    return memory


def save_results_to_csv(results, command, backend, binsize, protocol, had_errors=False):
    """Save benchmark results to a CSV file."""
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    if protocol:
        output_file = f"{output_dir}/{command}{backend}_{binsize}_{protocol}.csv"
    else:
        output_file = f"{output_dir}/{command}{backend}_{binsize}.csv"

    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Add error column if errors were detected
    if had_errors and "error" not in df.columns:
        df["error"] = "Error detected in file"

    # Save with NA representation for None values
    df.to_csv(output_file, index=False, na_rep="NA")

    print(f"Saved{'(with errors)' if had_errors else ''} to {output_file}")


def read_benchmark(file_path):
    print(f"Reading {file_path}")

    try:
        command, backend, binsize, protocol = extract_metadata_from_path(file_path)
    except ValueError as e:
        print(f"Error: {e}")
        return None

    results, had_errors = parse_benchmark_file(file_path)

    # Handle macOS memory regardless of errors
    try:
        if (
            "memory" in results
            and results["memory"]
            and not all(m is None for m in results["memory"])
        ):
            results["memory"] = handle_macos_memory(
                command, backend, protocol, results["memory"]
            )
    except AssertionError as e:
        print(f"Warning: Could not handle macOS memory: {e}")

    # Save results to CSV even if there were errors
    save_results_to_csv(results, command, backend, binsize, protocol, had_errors)

    return results


def parse_memory_from_logs(log_files):
    memory_values = []
    for log_file in log_files:
        with open(log_file, "r") as file:
            for line in file:
                if "maximum resident set size" in line:
                    memory_value = int(line.split()[0])
                    memory_values.append(memory_value)
                    break  # Go to next log file
    return memory_values


def apply_boxplot_style(bp, colors=["lightblue", "lightgreen"]):
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    plt.setp(bp["medians"], color="navy")
    plt.setp(bp["whiskers"], color="navy")
    plt.setp(bp["caps"], color="navy")
    plt.setp(bp["fliers"], marker="o", markerfacecolor="red", alpha=0.5)


def create_boxplot(
    ax,
    data1,
    data2,
    title=None,
    ylabel=None,
    is_first=True,
    labels=["Legacy", "Maturin"],
):
    """Helper function to create a boxplot on a given axis."""
    # Filter out None values
    valid_data1 = [x for x in data1 if x is not None]
    valid_data2 = [x for x in data2 if x is not None]

    # Handle case where there's no valid data
    if not valid_data1 and not valid_data2:
        ax.text(
            0.5,
            0.5,
            "No valid data available",
            ha="center",
            va="center",
            fontsize=14,
            color="red",
        )
        if is_first and ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title)
        return None

    # If one dataset is empty, use zeros and add text
    if not valid_data1:
        valid_data1 = [0]
        ax.text(
            0.25,
            0.05,
            "No valid data",
            ha="center",
            va="bottom",
            fontsize=10,
            color="red",
            transform=ax.transAxes,
        )

    if not valid_data2:
        valid_data2 = [0]
        ax.text(
            0.75,
            0.05,
            "No valid data",
            ha="center",
            va="bottom",
            fontsize=10,
            color="red",
            transform=ax.transAxes,
        )

    # Create boxplot
    bp = ax.boxplot([valid_data1, valid_data2], patch_artist=True, labels=labels)

    # Style the boxplot
    apply_boxplot_style(bp)

    # Add info about filtered values if needed
    orig_len1 = len(data1) if data1 else 0
    orig_len2 = len(data2) if data2 else 0
    valid_len1 = len(valid_data1) if valid_data1 != [0] else 0
    valid_len2 = len(valid_data2) if valid_data2 != [0] else 0

    if orig_len1 > valid_len1 or orig_len2 > valid_len2:
        ax.text(
            0.5,
            1.05,
            f"{valid_len1}/{orig_len1} vs {valid_len2}/{orig_len2} points",
            transform=ax.transAxes,
            ha="center",
            fontsize=8,
            color="blue",
        )

    if is_first and ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return bp


def make_boxplot(data1, data2, title, ylabel, labels=["Legacy", "Maturin"]):
    """Create a boxplot comparing two datasets."""
    fig, ax = plt.subplots(figsize=(8, 6))

    create_boxplot(ax, data1, data2, title, ylabel, labels=labels)

    plt.tight_layout()
    return fig


def make_protocol_boxplots(
    data1_files,
    data2_files,
    title,
    ylabel,
    metric="times",
    protocols=["ChIP", "RNA", "WGS"],
):
    """Create boxplots for multiple protocols."""
    fig, axes = plt.subplots(1, len(protocols), figsize=(15, 5))
    fig.suptitle(title)

    data1_file_list = data1_files.split(",")
    data2_file_list = data2_files.split(",")

    for idx, protocol in enumerate(protocols):
        # Get the appropriate axis
        if len(protocols) > 1:
            ax = axes[idx]
        else:
            ax = axes

        # Ensure we have enough files
        if idx >= len(data1_file_list) or idx >= len(data2_file_list):
            print(f"Warning: Not enough files for protocol {protocol}")
            ax.text(
                0.5,
                0.5,
                f"Missing data files for {protocol}",
                ha="center",
                va="center",
                fontsize=12,
                color="red",
            )
            ax.set_title(f"{protocol}")
            continue

        # Read the data
        data1 = read_benchmark(data1_file_list[idx])
        data2 = read_benchmark(data2_file_list[idx])

        # Extract the data for this metric, handling possible None values
        data1_metric = data1.get(metric, []) if data1 else []
        data2_metric = data2.get(metric, []) if data2 else []

        # Create the boxplot
        create_boxplot(
            ax,
            data1_metric,
            data2_metric,
            title=protocol,
            ylabel=ylabel if idx == 0 else None,
            is_first=(idx == 0),
        )

    plt.tight_layout()
    return fig


def process_multiprotocol(
    data1_files, data2_files, base_path, protocols, extension=".png"
):
    """Process multiple protocol files and generate comparison plots."""
    time_fig = make_protocol_boxplots(
        data1_files,
        data2_files,
        "Execution Time Comparison",
        "Time (s)",
        metric="times",
        protocols=protocols,
    )
    time_fig.savefig(f"{base_path}_time{extension}")

    mem_fig = make_protocol_boxplots(
        data1_files,
        data2_files,
        "Memory Usage Comparison",
        "Memory (MB)",
        metric="memory",
        protocols=protocols,
    )
    mem_fig.savefig(f"{base_path}_mem{extension}")


def process_single_files(data1_file, data2_file, base_path, extension=".png"):
    """Process single files and generate comparison plots."""
    data1 = read_benchmark(data1_file)
    data2 = read_benchmark(data2_file)

    # Handle different scenarios based on available data
    if data1 is None and data2 is None:
        print("Error: Both input files contain errors, cannot generate plots")
        return

    time_title = "Execution Time Comparison"
    mem_title = "Memory Usage Comparison"

    if data1 is None:
        print(f"Warning: Using only data from {data2_file}")
        data1_times = [None]
        data1_memory = [None]
        time_title += " (Partial Data)"
        mem_title += " (Partial Data)"
    else:
        data1_times = data1["times"]
        data1_memory = data1["memory"]

    if data2 is None:
        print(f"Warning: Using only data from {data1_file}")
        data2_times = [None]
        data2_memory = [None]
        time_title += " (Partial Data)"
        mem_title += " (Partial Data)"
    else:
        data2_times = data2["times"]
        data2_memory = data2["memory"]

    # Generate plots with available data
    time_fig = make_boxplot(data1_times, data2_times, time_title, "Time (s)")
    time_fig.savefig(f"{base_path}_time{extension}")

    mem_fig = make_boxplot(data1_memory, data2_memory, mem_title, "Memory (MB)")
    mem_fig.savefig(f"{base_path}_mem{extension}")


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="Generate benchmark comparison plots")

    parser.add_argument(
        "output_template", help="Output file template (e.g., 'output.png')"
    )
    parser.add_argument(
        "data1_files", help="First dataset file(s), comma-separated for protocol plots"
    )
    parser.add_argument(
        "data2_files", help="Second dataset file(s), comma-separated for protocol plots"
    )
    parser.add_argument(
        "--protocols", default="ChIP,RNA,WGS", help="Protocols for multi-protocol plots"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_command_line_args()
    protocols = args.protocols.split(",")
    output_template = args.output_template

    output_path = Path(output_template)
    base_path = output_path.with_suffix("")

    if "," in args.data1_files:
        process_multiprotocol(args.data1_files, args.data2_files, base_path, protocols)
    else:
        process_single_files(args.data1_files, args.data2_files, base_path)
