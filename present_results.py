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

    with open(file_path, "r") as file:
        _ = next(file)  # Skip header
        for line in file:
            parts = line.strip().split("\t")
            times.append(float(parts[0]))
            memory.append(float(parts[2]))
            max_uss.append(float(parts[4]))
            max_pss.append(float(parts[5]))
            io_in.append(float(parts[6]))
            io_out.append(float(parts[7]))
            mean_load.append(float(parts[8]))
            cpu_time.append(float(parts[9]))
            cpu_usage.append(float(parts[15].split()[0]))

    return {
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


def save_results_to_csv(results, command, backend, binsize, protocol):
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    if protocol:
        output_file = f"{output_dir}/{command}{backend}_{binsize}_{protocol}.csv"
    else:
        output_file = f"{output_dir}/{command}{backend}_{binsize}.csv"

    pd.DataFrame(results).to_csv(output_file, index=False)


def read_benchmark(file_path):
    print(f"Reading {file_path}")

    command, backend, binsize, protocol = extract_metadata_from_path(file_path)
    results = parse_benchmark_file(file_path)
    results["memory"] = handle_macos_memory(
        command, backend, protocol, results["memory"]
    )
    save_results_to_csv(results, command, backend, binsize, protocol)

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


def make_boxplot(data1, data2, title, ylabel):
    fig, ax = plt.subplots(figsize=(8, 5))
    bp = ax.boxplot([data1, data2], patch_artist=True, labels=["Legacy", "Maturin"])
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    apply_boxplot_style(bp)
    return fig


def make_protocol_boxplots(
    data1_files,
    data2_files,
    title,
    ylabel,
    metric="times",
    protocols=["ChIP", "RNA", "WGS"],
):
    fig, axes = plt.subplots(1, len(protocols), figsize=(15, 5))
    fig.suptitle(title)

    data1_file_list = data1_files.split(",")
    data2_file_list = data2_files.split(",")

    for idx, protocol in enumerate(protocols):
        # Ensure we have enough files
        if idx >= len(data1_file_list) or idx >= len(data2_file_list):
            print(f"Warning: Not enough files for protocol {protocol}")
            continue

        data1 = read_benchmark(data1_file_list[idx])
        data2 = read_benchmark(data2_file_list[idx])

        bp = axes[idx].boxplot(
            [data1[metric], data2[metric]],
            patch_artist=True,
            labels=["Legacy", "Maturin"],
        )
        axes[idx].set_title(f"{protocol}")
        axes[idx].set_ylabel(ylabel if idx == 0 else "")
        apply_boxplot_style(bp)

    plt.tight_layout()
    return fig


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

    # Fix file extension handling with Path
    output_path = Path(output_template)
    base_path = output_path.with_suffix("")
    extension = output_path.suffix

    if "," in args.data1_files:
        # Multi-protocol comparison
        time_fig = make_protocol_boxplots(
            args.data1_files,
            args.data2_files,
            "Execution Time Comparison",
            "Time (s)",
            metric="times",
            protocols=protocols,
        )
        time_fig.savefig(f"{base_path}_time{extension}")

        mem_fig = make_protocol_boxplots(
            args.data1_files,
            args.data2_files,
            "Memory Usage Comparison",
            "Memory (MB)",
            metric="memory",
            protocols=protocols,
        )
        mem_fig.savefig(f"{base_path}_mem{extension}")
    else:
        # Single file comparison
        data1 = read_benchmark(args.data1_files)
        data2 = read_benchmark(args.data2_files)

        time_fig = make_boxplot(
            data1["times"], data2["times"], "Execution Time Comparison", "Time (s)"
        )
        time_fig.savefig(f"{base_path}_time{extension}")

        mem_fig = make_boxplot(
            data1["memory"], data2["memory"], "Memory Usage Comparison", "Memory (MB)"
        )
        mem_fig.savefig(f"{base_path}_mem{extension}")
