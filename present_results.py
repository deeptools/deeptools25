#!/usr/bin/env python
import argparse
import datetime
import glob
import logging
import os
import platform
import re
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Global constants for plot labels
PREV_VERSION_LABEL = "Previous version"
CURRENT_VERSION_LABEL = "4.0"


# Matplotlib settings for our use-case(s)
mpl.use("Agg")

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from matplotlib.cbook import MatplotlibDeprecationWarning

    warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)
except ImportError:
    warnings.filterwarnings("ignore", message=".*deprecated.*", category=Warning)


# Set up logging
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = f"logs/present_results_{timestamp}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(),  # Console output
    ],
)

logger = logging.getLogger("benchmark")
logger.info(f"Logging to file: {log_file}")

# Global variable to track the first command name encountered
cmd_being_process = None


class BenchmarkResult:
    """
    Represents benchmark data with standardized error handling.
    Separates data from Python (Snakemake) and kernel (time command).
    """

    def __init__(
        self, python_data=None, kernel_data=None, error=None, error_indices=None
    ):
        """
        Initialize benchmark result with data from different sources and optional error.

        Args:
            python_data: Dictionary of metrics from Python/Snakemake
            kernel_data: Dictionary of metrics from time command
            error: Error information (string or Exception)
            error_indices: List of indices where errors were detected
        """
        self.python_data = python_data or {}
        self.kernel_data = kernel_data or {}
        self.error = error
        self.error_indices = error_indices or []

    @property
    def data(self):
        """
        Combined data dictionary for backward compatibility.
        Prioritizes kernel data over Python data when both exist.
        """
        combined = {}
        combined.update(self.python_data)
        combined.update(self.kernel_data)
        return combined

    @property
    def is_valid(self):
        """Check if result contains valid data despite possible errors."""
        return any(
            bool(data and any(len(v) > 0 for v in data.values() if isinstance(v, list)))
            for data in [self.python_data, self.kernel_data]
        )

    @property
    def is_error(self):
        """Check if result has an error."""
        return self.error is not None

    @property
    def is_empty(self):
        """Check if result has no usable data."""
        return not self.is_valid

    def get_metric(self, metric_name, source="combined"):
        """
        Safely get a metric value list from the specified data source.

        Args:
            metric_name: Name of the metric to retrieve
            source: Data source - 'python', 'kernel', or 'combined' (default)

        Returns:
            List of metric values or empty list if not found
        """
        if source == "python":
            value = self.python_data.get(metric_name, [])
        elif source == "kernel":
            value = self.kernel_data.get(metric_name, [])
        else:  # combined - prefer kernel data when available
            value = self.kernel_data.get(
                metric_name, self.python_data.get(metric_name, [])
            )

        return value if isinstance(value, list) else []

    def mark_index_as_error(self, index):
        """
        Mark a specific measurement index as invalid in all metrics.
        This helps when we detect an error in one measurement but want to
        invalidate that measurement across all metrics.

        Args:
            index: Index to mark as invalid (set to None)
        """
        if index in self.error_indices:
            return  # Already marked

        self.error_indices.append(index)

        # Mark as None in Python data
        for metric, values in self.python_data.items():
            if isinstance(values, list) and 0 <= index < len(values):
                values[index] = None

        # Mark as None in kernel data
        for metric, values in self.kernel_data.items():
            if isinstance(values, list) and 0 <= index < len(values):
                values[index] = None


def validate_measurements(results, expected_count):
    """
    Validate that we have the expected number of measurements or report discrepancies.
    Checks both Python and kernel data sources independently.

    Args:
        results: BenchmarkResult object
        expected_count: Expected number of measurements

    Returns:
        Updated BenchmarkResult object with validation info
    """
    if not results or not results.is_valid or expected_count is None:
        return results

    # Check counts for each data source
    validation_issues = {}

    # Validate Python data first
    for metric, values in results.python_data.items():
        actual_count = len(values) if values is not None else 0
        if actual_count != expected_count:
            validation_issues[f"python_{metric}"] = (
                f"Expected {expected_count} values, got {actual_count}"
            )

    # Then validate kernel data
    for metric, values in results.kernel_data.items():
        actual_count = len(values) if values is not None else 0
        if actual_count != expected_count:
            validation_issues[f"kernel_{metric}"] = (
                f"Expected {expected_count} values, got {actual_count}"
            )

    # Add validation results to the result object without marking as error
    if validation_issues:
        results.data["validation_issues"] = validation_issues
        logger.warning(f"Validation issues detected: {validation_issues}")

        # Only set as error if there are serious inconsistencies (missing measurements)
        critical_issues = any(
            "got 0" in str(issue) for issue in validation_issues.values()
        )
        if critical_issues and not results.is_error:
            results.error = f"Missing measurements detected: {validation_issues}"

    return results


def extract_metadata_from_path(file_path):
    try:
        match = re.search(
            r"/metrics_(\w+?)(\d+)_(?:.*?)(?:_bs(\d+))?(?:_(\w+))?\.txt", file_path
        )
        command = match.group(1)  # 'bamCompare' or 'bamCoverage'
        backend = match.group(2)  # '1'
        binsize = match.group(3)  # '100'
        protocol = match.group(4)  # 'chip' or None
        return command, backend, binsize, protocol
    except (AttributeError, IndexError) as e:
        raise ValueError(f"Could not parse filename pattern from {file_path}: {e}")


def parse_benchmark_file(file_path, n_threads, platform_override="auto"):
    """
    Parse benchmark files, separating Snakemake benchmark (python_data) and time command output (kernel_data).

    Wall time is the only metric that would come primarily from Snakemake.

    Args:
        file_path: Path to benchmark log file
        n_threads: Number of threads used
        platform_override: Override platform detection (auto, Linux, Darwin)

    Returns:
        BenchmarkResult with separated python and kernel data
    """
    # Get platform with possible override
    current_platform = get_platform(platform_override)

    # Initialize data arrays
    python_wall_times, python_max_rss = [], []
    error_detected = False
    error_message = None
    error_indices = []

    # First - parse Snakemake benchmark file (metrics_*.txt)
    try:
        with open(file_path, "r") as file:
            content = file.readlines()

            # Check for global errors in the entire content
            full_content = "".join(content)

            # Enhanced error detection - check for various error patterns
            error_patterns = {
                "does not exist": "File not found error",
                "Read block operation failed": "Read block operation failed",
                "BamTruncatedRecord": "BAM file truncated",
                "Error parsing record": "Error parsing BAM record",
                "Command exited with non-zero status": "Command failed",
                "segmentation fault": "Segmentation fault",
                "Killed": "Process was killed",
                "out of memory": "Out of memory error",
                "cannot allocate memory": "Memory allocation failed",
                "Traceback": "Python exception occurred",
                "panicked at": "Rust panic occurred",
                "Exception": "Exception occurred",
                "error:": "Error detected",
                "failed with exit code": "Process failed",
                "No such file or directory": "File not found",
                "Permission denied": "Permission denied",
            }

            # Identify which run number this error applies to
            current_run_index = -1

            # Look for iteration markers to determine run number
            for line in content:
                if "iter_" in line and ".bw" in line:
                    try:
                        # Try to extract run number from command line with iter_X.bw
                        run_match = re.search(r"iter_(\d+)\.bw", line)
                        if run_match:
                            run_num = int(run_match.group(1))
                            current_run_index = run_num - 1  # 0-based index
                    except (ValueError, IndexError):
                        pass

            # Check for each known error pattern
            for pattern, msg in error_patterns.items():
                if pattern in full_content.lower():
                    error_detected = True
                    error_message = msg
                    logger.warning(f"{msg} found in {file_path}")

                    # Mark the specific run as having an error if we know which one
                    if (
                        current_run_index >= 0
                        and current_run_index not in error_indices
                    ):
                        error_indices.append(current_run_index)

            # Skip header line for Snakemake benchmark data
            content = content[1:]

            # Parse Snakemake benchmark data
            run_index = 0
            for line_num, line in enumerate(content, 2):
                # Skip error message lines
                if any(
                    skip in line
                    for skip in ["does not exist", "Error", "Traceback", "panicked"]
                ):
                    continue

                try:
                    parts = line.strip().split("\t")

                    # Skip lines with incomplete data
                    if len(parts) < 3:
                        continue

                    # Parse wall time (column 1) and max_rss (column 3)
                    wall_time = float(parts[0])
                    max_rss_val = float(parts[2])

                    # If this run had an error, mark the data as invalid
                    if run_index in error_indices:
                        wall_time = None
                        max_rss_val = None

                    # Store python data
                    python_wall_times.append(wall_time)
                    python_max_rss.append(max_rss_val)

                    # Move to next run
                    run_index += 1

                except (IndexError, ValueError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

    except FileNotFoundError:
        error_detected = True
        error_message = "File not found"
        logger.warning(f"{error_message}: {file_path}")

    # Now parse time command output from log files
    kernel_data = {
        "cpu_time": [],
        "wall_time": [],
        "max_rss": [],
        "max_uss": [],
        "max_pss": [],
        "io_in": [],
        "io_out": [],
        "mean_load": [],
        "cpu_usage_pct": [],
    }

    # Get command name, backend, binsize, protocol from file path
    try:
        command, backend, binsize, protocol = extract_metadata_from_path(file_path)
    except Exception as e:
        logger.warning(f"Could not extract metadata from file path: {e}")
        command = "unknown"
        backend = "unknown"
        protocol = None

    # The log files follow a different naming pattern than the benchmark files
    # Format should be: logs/bamCompare1_1.txt, logs/bamCoverage1_chip_1.txt etc.
    log_files = []

    # Try different log patterns based on command structure
    if protocol:
        # For commands like bamCoverage that have a protocol
        log_pattern = f"logs/{command}{backend}_{protocol}_[0-9]*.txt"
    else:
        # For commands like bamCompare that don't have a protocol
        log_pattern = f"logs/{command}{backend}_[0-9]*.txt"

    log_files = glob.glob(log_pattern)

    if not log_files:
        logger.warning(
            f"No log files found for {command}{backend} using pattern: {log_pattern}"
        )

    # Log what we found for debugging
    logger.debug(
        f"Found {len(log_files)} log files for {command}{backend} at: {log_files}"
    )

    # Now instead of using platform.system(), use the current_platform variable:
    if current_platform == "Linux":
        # Extract time command data from Linux output
        for log_file in log_files:
            try:
                with open(log_file, "r") as f:
                    log_content = f.read()

                    # Extract run number from log file name
                    run_match = re.search(r"_(\d+)\.txt$", log_file)
                    run_index = int(run_match.group(1)) - 1 if run_match else -1

                    # Skip if this run had errors
                    if run_index in error_indices:
                        kernel_data["cpu_time"].append(None)
                        kernel_data["wall_time"].append(None)
                        kernel_data["max_rss"].append(None)
                        kernel_data["max_uss"].append(None)
                        kernel_data["max_pss"].append(None)
                        kernel_data["io_in"].append(None)
                        kernel_data["io_out"].append(None)
                        kernel_data["mean_load"].append(None)
                        kernel_data["cpu_usage_pct"].append(None)
                        continue

                    # Parse Linux time command output
                    user_match = re.search(
                        r"User time \(seconds\): (\d+\.\d+)", log_content
                    )
                    sys_match = re.search(
                        r"System time \(seconds\): (\d+\.\d+)", log_content
                    )
                    wall_match = re.search(
                        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\d+:\d+\.\d+)",
                        log_content,
                    )
                    memory_match = re.search(
                        r"Maximum resident set size \(kbytes\): (\d+)", log_content
                    )

                    # Get CPU time (user + system)
                    if user_match and sys_match:
                        user_time = float(user_match.group(1))
                        sys_time = float(sys_match.group(1))
                        cpu_time = user_time + sys_time
                    else:
                        cpu_time = None

                    if wall_match:
                        # Convert mm:ss to seconds
                        time_str = wall_match.group(1)
                        minutes, seconds = time_str.split(":")
                        wall_time = float(minutes) * 60 + float(seconds)
                    else:
                        wall_time = None

                    if memory_match:
                        # Convert to MB
                        memory = float(memory_match.group(1)) / 1024
                    else:
                        memory = None

                    # Append kernel data
                    kernel_data["cpu_time"].append(cpu_time)
                    kernel_data["wall_time"].append(wall_time)
                    kernel_data["max_rss"].append(memory)
                    # Use None for metrics we can't parse directly from time output
                    kernel_data["max_uss"].append(None)
                    kernel_data["max_pss"].append(None)
                    kernel_data["io_in"].append(None)
                    kernel_data["io_out"].append(None)
                    kernel_data["mean_load"].append(None)
                    kernel_data["cpu_usage_pct"].append(None)

            except Exception as e:
                logger.warning(f"Error parsing kernel metrics from {log_file}: {e}")

    elif current_platform == "Darwin":
        # Extract time command data from macOS output
        for log_file in log_files:
            try:
                with open(log_file, "r") as f:
                    log_content = f.read()

                    # Extract run number from log file name
                    run_match = re.search(r"_(\d+)\.txt$", log_file)
                    run_index = int(run_match.group(1)) - 1 if run_match else -1

                    # Skip if this run had errors
                    if run_index in error_indices:
                        kernel_data["cpu_time"].append(None)
                        kernel_data["wall_time"].append(None)
                        kernel_data["max_rss"].append(None)
                        kernel_data["max_uss"].append(None)
                        kernel_data["max_pss"].append(None)
                        kernel_data["io_in"].append(None)
                        kernel_data["io_out"].append(None)
                        kernel_data["mean_load"].append(None)
                        kernel_data["cpu_usage_pct"].append(None)
                        continue

                    # Parse macOS time command output
                    real_match = re.search(r"real\s+(\d+\.\d+)", log_content)
                    user_match = re.search(r"user\s+(\d+\.\d+)", log_content)
                    sys_match = re.search(r"sys\s+(\d+\.\d+)", log_content)
                    memory_match = re.search(r"(\d+)\s+maximum resident", log_content)

                    if user_match and sys_match:
                        user_time = float(user_match.group(1))
                        sys_time = float(sys_match.group(1))
                        cpu_time = user_time + sys_time
                    else:
                        cpu_time = None

                    if real_match:
                        wall_time = float(real_match.group(1))
                    else:
                        wall_time = None

                    if memory_match:
                        # Convert to MB
                        memory = float(memory_match.group(1)) / 1024
                    else:
                        memory = None

                    # Append kernel data
                    kernel_data["cpu_time"].append(cpu_time)
                    kernel_data["wall_time"].append(wall_time)
                    kernel_data["max_rss"].append(memory)
                    # Use None for metrics we can't parse directly from time output
                    kernel_data["max_uss"].append(None)
                    kernel_data["max_pss"].append(None)
                    kernel_data["io_in"].append(None)
                    kernel_data["io_out"].append(None)
                    kernel_data["mean_load"].append(None)
                    kernel_data["cpu_usage_pct"].append(None)

            except Exception as e:
                logger.warning(f"Error parsing kernel metrics from {log_file}: {e}")

    # Create Python data dictionary with clearer naming
    python_data = {
        "wall_time": python_wall_times,
        "memory_mb": python_max_rss,
    }

    # Calculate derived metrics if we have enough data
    if (
        n_threads
        and any(x is not None for x in kernel_data["cpu_time"])
        and any(x is not None for x in kernel_data["wall_time"])
    ):
        # Calculate parallelization efficiency
        parallelization_efficiency = calculate_parallelization_efficiency(
            kernel_data["wall_time"], kernel_data["cpu_time"], n_threads
        )
        kernel_data["parallelization_efficiency"] = parallelization_efficiency

    # Add memory metrics if available
    if any(x is not None for x in kernel_data["max_rss"]) and any(
        x is not None for x in kernel_data["max_uss"]
    ):
        memory_metrics = calculate_memory_metrics(
            kernel_data["max_rss"], kernel_data["max_uss"], n_threads
        )
        kernel_data.update(memory_metrics)

    # If we detected errors and got no data, add placeholder values
    if error_detected and not python_wall_times:
        python_data = {key: [None] for key in python_data}
        kernel_data = {key: [None] for key in kernel_data}

    return BenchmarkResult(
        python_data=python_data,
        kernel_data=kernel_data,
        error=error_message if error_detected else None,
        error_indices=error_indices,
    )


def handle_macos_memory(command, backend, protocol, memory, platform_override="auto"):
    # Use the overridden platform
    if get_platform(platform_override) != "Darwin":
        return memory

    n_times = len(memory)

    try:
        if sum(x or 0 for x in memory) != 0:  # Allow for None values
            return memory
    except TypeError:
        logger.warning("Memory values contain non-numeric values")
        return memory

    try:
        if command == "bamCoverage":
            if protocol is None:
                logger.warning(
                    f"Couldn't retrieve protocol for {command}{backend}, skipping macOS memory handling"
                )
                return memory
            log_files = glob.glob(
                f"logs/metrics_{command}{backend}_{protocol}_[0-9]*.txt"
            )
        else:
            log_files = glob.glob(f"logs/metrics_{command}{backend}_[0-9]*.txt")

        if not log_files:
            logger.warning(
                f"No log files found for macOS memory handling: {command}{backend}"
            )
            return memory

        memory = parse_memory_from_logs(log_files)

        if len(memory) != n_times:
            logger.warning(f"Expected {n_times} memory values but got {len(memory)}")
            # Pad with None if needed
            if len(memory) < n_times:
                memory.extend([None] * (n_times - len(memory)))
            # Truncate if too many
            elif len(memory) > n_times:
                memory = memory[:n_times]

        return memory
    except Exception as e:
        logger.warning(f"Error handling macOS memory: {e}")
        return memory


def save_results_to_csv(results, command, backend, binsize, protocol):
    """
    Save benchmark results to a CSV file with clear separation of Python and kernel metrics.

    Args:
        results: BenchmarkResult object containing the data
        command: Tool command name (e.g., 'bamCoverage')
        backend: Backend identifier (e.g., '1' or '2')
        binsize: Bin size used
        protocol: Protocol name or None
    """
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate file name
    if protocol:
        output_file = f"{output_dir}/{command}{backend}_{binsize}_{protocol}.csv"
    else:
        output_file = f"{output_dir}/{command}{backend}_{binsize}.csv"

    # Create a dataframe with combined and aligned Python and kernel metrics
    data = {}

    # Get number of data points
    max_len = 0
    for source, data_dict in [
        ("python", results.python_data),
        ("kernel", results.kernel_data),
    ]:
        for key, values in data_dict.items():
            if isinstance(values, list):
                max_len = max(max_len, len(values))

    # Create run index
    data["run_index"] = list(range(1, max_len + 1))

    # Add Python (Snakemake) metrics with clear prefixes
    for key, values in results.python_data.items():
        if isinstance(values, list):
            # Ensure all lists have the same length
            padded_values = values + [None] * (max_len - len(values))
            column_name = f"python_{key}"
            data[column_name] = padded_values

    # Add kernel (time command) metrics with clear prefixes
    for key, values in results.kernel_data.items():
        if isinstance(values, list):
            # Ensure all lists have the same length
            padded_values = values + [None] * (max_len - len(values))
            column_name = f"kernel_{key}"
            data[column_name] = padded_values

    # Add metadata columns
    data["command"] = [command] * max_len
    data["backend"] = [backend] * max_len
    data["bin_size"] = [binsize] * max_len
    if protocol:
        data["protocol"] = [protocol] * max_len

    # Add error info if present
    if results.is_error:
        data["error_message"] = [results.error] * max_len

    # Add error indices if any
    if hasattr(results, "error_indices") and results.error_indices:
        # Mark which runs had errors
        data["has_error"] = [
            "Yes" if i in results.error_indices else "No" for i in range(max_len)
        ]

    # Create and save DataFrame
    df = pd.DataFrame(data)

    # Set columns in a logical order
    cols = ["run_index"]
    if "has_error" in df.columns:
        cols.append("has_error")

    # Time metrics (wall time first, then CPU time)
    time_metrics = [col for col in df.columns if "time" in col]
    wall_time_metrics = [col for col in time_metrics if "wall" in col]
    cpu_time_metrics = [col for col in time_metrics if "cpu" in col]
    other_time_metrics = [
        col for col in time_metrics if col not in wall_time_metrics + cpu_time_metrics
    ]
    cols.extend(wall_time_metrics + cpu_time_metrics + other_time_metrics)

    # Memory metrics
    memory_metrics = [
        col
        for col in df.columns
        if "memory" in col
        or "mem" in col
        or "rss" in col
        or "uss" in col
        or "pss" in col
    ]
    cols.extend(memory_metrics)

    # Other metrics
    other_metrics = [
        col
        for col in df.columns
        if col not in cols
        and col not in ["command", "backend", "bin_size", "protocol", "error_message"]
    ]
    cols.extend(other_metrics)

    # Metadata columns at the end
    meta_cols = [
        col
        for col in ["command", "backend", "bin_size", "protocol", "error_message"]
        if col in df.columns
    ]
    cols.extend(meta_cols)

    # Reorder columns
    df = df[cols]

    # Add some calculated columns for easy comparison
    if "python_wall_time" in df.columns and "kernel_cpu_time" in df.columns:
        non_null_mask = (
            df["python_wall_time"].notnull() & df["kernel_cpu_time"].notnull()
        )
        if non_null_mask.any():
            # CPU/Wall ratio
            df.loc[non_null_mask, "cpu_vs_wall_ratio"] = (
                df.loc[non_null_mask, "kernel_cpu_time"]
                / df.loc[non_null_mask, "python_wall_time"]
            )

            # Add workload classification if we know thread count
            n_threads = results.data.get("n_threads")
            if n_threads:
                df.loc[non_null_mask, "expected_ratio"] = n_threads
                df.loc[non_null_mask, "ratio_to_expected"] = (
                    df.loc[non_null_mask, "cpu_vs_wall_ratio"] / n_threads
                )

                # Add workload classification
                def classify_workload(ratio):
                    if ratio < 0.5:
                        return "I/O-bound"
                    elif 0.7 <= ratio <= 1.3:
                        return "CPU-bound"
                    elif ratio > 1.3:
                        return "Super-linear"
                    else:
                        return "Mixed"

                if "ratio_to_expected" in df.columns:
                    df.loc[non_null_mask, "workload_type"] = df.loc[
                        non_null_mask, "ratio_to_expected"
                    ].apply(classify_workload)

    # Save with NA representation for None values
    df.to_csv(output_file, index=False, na_rep="NA")

    logger.info(f"Saved{'(with errors)' if results.is_error else ''} to {output_file}")

    return output_file


def read_benchmark(
    file_path, expected_count=None, n_threads=None, platform_override="auto"
):
    """
    Read and parse a benchmark file with validation against expected measurement count.

    Args:
        file_path: Path to the benchmark file
        expected_count: Expected number of measurements
        n_threads: Number of threads used in benchmark
        platform_override: Override platform detection (auto, Linux, Darwin)

    Returns:
        BenchmarkResult object
    """
    global cmd_being_process  # Use the global variable
    logger.info(f"Reading {file_path}")

    try:
        command, backend, binsize, protocol = extract_metadata_from_path(file_path)

        if cmd_being_process is None:
            cmd_being_process = command

    except ValueError as e:
        logger.error(f"{e}")
        return BenchmarkResult(error=str(e))

    result = parse_benchmark_file(file_path, n_threads, platform_override)

    # Handle macOS memory - use overridden platform detection here too
    try:
        memory_data = result.get_metric("memory")
        if memory_data and not all(m is None for m in memory_data):
            updated_memory = handle_macos_memory(
                command,
                backend,
                protocol,
                memory_data,
                platform_override=platform_override,
            )
            result.data["memory"] = updated_memory
    except Exception as e:
        logger.warning(f"Could not handle macOS memory: {e}")

    # Validate against expected measurement count if provided
    if expected_count is not None:
        result = validate_measurements(result, expected_count)

    # Check consistency between Python and kernel data (only for Linux, unless overridden)
    here_platform = get_platform(platform_override)
    if here_platform == "Linux":
        inconsistencies = validate_data_consistency(
            result, n_threads, current_platform=here_platform
        )
        if inconsistencies:
            # Only log a warning if we have actual inconsistencies (outliers in CPU-bound tasks)
            if any(key.startswith("cpu_time_") for key in inconsistencies):
                # Only log a summary instead of each individual measurement
                inconsistent_measurements = len(
                    [k for k in inconsistencies if k.startswith("cpu_time_")]
                )
                total_measurements = len(result.get_metric("wall_time"))
                logger.warning(
                    f"CPU/Wall time inconsistencies in {inconsistent_measurements}/{total_measurements} "
                    f"measurements in {file_path}"
                )
                # Log just the first example as a sample
                first_key = next(
                    k for k in inconsistencies if k.startswith("cpu_time_")
                )
                logger.warning(f"Example: {first_key}: {inconsistencies[first_key]}")
            elif "summary" in inconsistencies:
                # For non-CPU-bound tasks with significant deviation, log as warning
                logger.warning(
                    f"Workload anomaly detected: {inconsistencies['summary']}"
                )

            # Add to result data but don't treat as errors
            result.data["data_inconsistencies"] = inconsistencies

    # Save results to CSV even if there were errors
    save_results_to_csv(result, command, backend, binsize, protocol)

    return result


def parse_memory_from_logs(log_files):
    memory_values = []
    for log_file in log_files:
        try:
            with open(log_file, "r") as file:
                for line in file:
                    if "maximum resident set size" in line:
                        memory_value = int(line.split()[0])
                        memory_values.append(memory_value)
                        break  # Go to next log file
        except Exception as e:
            logger.warning(f"Error reading log file {log_file}: {e}")
    return memory_values


def validate_data_consistency(result, n_threads=None, current_platform=None):
    """
    Validate consistency between wall time and CPU time measurements.

    Args:
        result: BenchmarkResult object
        n_threads: Number of threads used (if known)
        current_platform: Platform detection passed from read_benchmark

    Returns:
        Dictionary of inconsistencies found
    """

    # Skip validation for non-Linux platforms
    assert current_platform is not None
    if current_platform != "Linux" or not n_threads or n_threads <= 1:
        return {}

    inconsistencies = {}

    # Get wall time and CPU time values
    wall_times = result.get_metric("wall_time", "python")
    cpu_times = result.get_metric("cpu_time", "kernel")

    # Calculate average ratio to detect patterns
    valid_ratios = []
    for wall, cpu in zip(wall_times, cpu_times):
        if wall is not None and cpu is not None and wall > 0:
            valid_ratios.append(cpu / wall)

    if not valid_ratios:
        return {}

    # Calculate median ratio and variance
    median_ratio = np.median(valid_ratios)

    # Store thread count in the result data for CSV output
    result.data["n_threads"] = n_threads

    # Determine if workload is CPU-bound (ratio close to thread count)
    ratio_of_expected = median_ratio / n_threads
    cpu_bound = 0.7 <= ratio_of_expected <= 1.3

    # Store workload analysis in result data
    workload_type = (
        "CPU-bound"
        if cpu_bound
        else "I/O-bound"
        if ratio_of_expected < 0.7
        else "Super-linear"
    )

    # Always log workload analysis as INFO
    logger.info(
        f"Workload analysis: Task appears to be {workload_type} - "
        f"median CPU/Wall ratio ({median_ratio:.1f}) is {abs(ratio_of_expected - 1.0) * 100:.0f}% "
        f"{'higher' if ratio_of_expected > 1 else 'lower'} than expected for {n_threads} threads"
    )

    # Store workload analysis in result
    result.data["workload_analysis"] = {
        "type": workload_type,
        "median_ratio": median_ratio,
        "expected_ratio": n_threads,
        "ratio_of_expected": ratio_of_expected,
    }

    # For CPU-bound tasks, flag individual outliers as INFO level observations
    if cpu_bound:
        for i, (wall, cpu) in enumerate(zip(wall_times, cpu_times)):
            if wall is None or cpu is None or wall == 0:
                continue

            ratio = cpu / wall

            # Only flag significant outliers - be more lenient for longer-running tasks
            threshold = 0.5 if wall < 10 else 0.4 if wall < 60 else 0.3

            if abs(ratio - median_ratio) / median_ratio > threshold:
                message = (
                    f"Measurement {i + 1}: CPU/Wall ratio ({ratio:.1f}) differs from median "
                    f"({median_ratio:.1f}) by {abs(ratio - median_ratio) / median_ratio * 100:.0f}%"
                )
                # Store in results but don't warn - just info
                inconsistencies[f"cpu_time_{i}"] = message

                # Log as INFO instead of WARNING
                logger.info(f"Outlier detected: {message}")

    # For non-CPU-bound tasks, use WARNING only for serious inconsistencies
    elif abs(ratio_of_expected - 1.0) > 1.5:  # Increased threshold to 150% for warnings
        description = "higher" if ratio_of_expected > 1 else "lower"

        # Create a single summary as a warning for severe anomalies
        warning_msg = (
            f"SIGNIFICANT ANOMALY: Task appears to be {workload_type} - "
            f"median CPU/Wall ratio ({median_ratio:.1f}) is {abs(ratio_of_expected - 1.0) * 100:.0f}% {description} "
            f"than expected for {n_threads} threads"
        )
        inconsistencies["summary"] = warning_msg
        logger.warning(warning_msg)

    # For moderate inconsistencies, use INFO level (100-150% different from expected)
    elif abs(ratio_of_expected - 1.0) > 1.0:
        description = "higher" if ratio_of_expected > 1 else "lower"

        # Create a single summary as an info message
        info_msg = (
            f"MODERATE ANOMALY: Task appears to be {workload_type} - "
            f"median CPU/Wall ratio ({median_ratio:.1f}) is {abs(ratio_of_expected - 1.0) * 100:.0f}% {description} "
            f"than expected for {n_threads} threads"
        )
        inconsistencies["summary"] = info_msg
        logger.info(info_msg)

    return inconsistencies


def calculate_parallelization_efficiency(wall_times, cpu_times, n_threads):
    """
    Calculate parallelization efficiency as (CPU time / Wall time) / n_threads

    Args:
        wall_times: List of wall times
        cpu_times: List of CPU times
        n_threads: Number of threads used

    Returns:
        List of parallelization efficiency values
    """
    efficiencies = []
    for wall, cpu in zip(wall_times, cpu_times):
        if wall is not None and cpu is not None and wall > 0 and n_threads > 0:
            efficiency = (cpu / wall) / n_threads
            efficiencies.append(efficiency)
        else:
            efficiencies.append(None)
    return efficiencies


def calculate_memory_metrics(max_rss, max_uss, n_threads):
    """
    Calculate advanced memory metrics

    Args:
        max_rss: List of maximum RSS values
        max_uss: List of maximum USS values
        n_threads: Number of threads used

    Returns:
        Dictionary containing memory efficiency metrics
    """
    results = {
        "memory_per_thread": [],  # Memory used per thread
        "uss_rss_ratio": [],  # How efficiently memory is shared
        "memory_efficiency": [],  # Approximation of memory scaling efficiency
    }

    for rss, uss in zip(max_rss, max_uss):
        if rss is not None and uss is not None and rss > 0:
            # Memory used per thread (RSS)
            mem_per_thread = rss / n_threads
            results["memory_per_thread"].append(mem_per_thread)

            # USS/RSS ratio (closer to 1 means less memory sharing)
            uss_rss = uss / rss
            results["uss_rss_ratio"].append(uss_rss)

            # Memory efficiency approximation
            # Higher values may indicate better memory usage patterns
            # This is a heuristic - memory scaling is complex
            uss_per_thread = uss / n_threads
            mem_efficiency = 1.0 - (uss_per_thread / rss)
            results["memory_efficiency"].append(mem_efficiency)
        else:
            results["memory_per_thread"].append(None)
            results["uss_rss_ratio"].append(None)
            results["memory_efficiency"].append(None)

    return results


def apply_boxplot_style(bp, colors=["lightblue", "lightgreen"]):
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    plt.setp(bp["medians"], color="navy")
    plt.setp(bp["whiskers"], color="navy")
    plt.setp(bp["caps"], color="navy")
    plt.setp(bp["fliers"], marker="o", markerfacecolor="red", alpha=0.5)


def process_single_files(
    data1_file,
    data2_file,
    base_path,
    n_times,
    n_threads,
    extension=".png",
    platform_override="auto",
):
    """Process single files and generate comparison plots with both CPU and wall time."""
    result1 = read_benchmark(data1_file, n_times, n_threads, platform_override)
    result2 = read_benchmark(data2_file, n_times, n_threads, platform_override)

    # Set up titles with error indicators if needed
    error_suffix = (
        " (Errors)"
        if (result1 and result1.is_error) or (result2 and result2.is_error)
        else ""
    )

    # Generate CPU time plot
    cpu_title = f"CPU Time Comparison{error_suffix}"

    # Handle different scenarios based on available data
    if result1 is None and result2 is None:
        logger.error("Both input files contain errors, cannot generate plots")
        return

    # Set up titles with error indicators if needed
    error_suffix = (
        " (Errors)"
        if (result1 and result1.is_error) or (result2 and result2.is_error)
        else ""
    )

    # Generate CPU time plot
    cpu_title = f"CPU Time Comparison{error_suffix}"
    cpu_time_fig = make_boxplot(
        result1, result2, cpu_title, "CPU Time (s)", metric="times", use_cpu_time=True
    )
    cpu_time_fig.savefig(f"{base_path}_cputime{extension}")
    logger.info(f"Saved CPU time plot to {base_path}_cputime{extension}")

    # Generate wall time plot
    wall_title = f"Wall Time Comparison{error_suffix}"
    wall_time_fig = make_boxplot(
        result1,
        result2,
        wall_title,
        "Wall Time (s)",
        metric="times",
        use_cpu_time=False,
    )
    wall_time_fig.savefig(f"{base_path}_walltime{extension}")
    logger.info(f"Saved wall time plot to {base_path}_walltime{extension}")

    # Generate memory plot (unchanged)
    mem_title = f"Memory Usage Comparison{error_suffix}"
    mem_fig = make_boxplot(result1, result2, mem_title, "Memory (MB)", metric="memory")
    mem_fig.savefig(f"{base_path}_memory{extension}")
    logger.info(f"Saved memory plot to {base_path}_memory{extension}")

    # Generate efficiency metrics plots
    if result1 and result1.is_valid and result2 and result2.is_valid:
        # Extract efficiency metrics
        efficiency_metrics1 = {
            "parallelization_efficiency": result1.get_metric(
                "parallelization_efficiency", "kernel"
            ),
            "memory_efficiency": result1.get_metric("memory_efficiency", "kernel"),
            "uss_rss_ratio": result1.get_metric("uss_rss_ratio", "kernel"),
        }

        efficiency_metrics2 = {
            "parallelization_efficiency": result2.get_metric(
                "parallelization_efficiency", "kernel"
            ),
            "memory_efficiency": result2.get_metric("memory_efficiency", "kernel"),
            "uss_rss_ratio": result2.get_metric("uss_rss_ratio", "kernel"),
        }

        # Add efficiency metrics plot
        plot_efficiency_metrics(
            efficiency_metrics1,
            efficiency_metrics2,
            title=f"Efficiency Metrics{error_suffix}",
            output_prefix=f"{base_path}_efficiency",
        )

        # Add workload analysis plot
        plot_workload_analysis(
            result1,
            result2,
            title=f"Workload Analysis{error_suffix}",
            output_prefix=base_path,
            n_threads=n_threads,
        )


def process_multiprotocol(
    data1_files,
    data2_files,
    base_path,
    protocols,
    n_times,
    n_threads,
    extension=".png",
    platform_override="auto",
):
    """Process multiple protocol files and generate comparison plots with both CPU and wall time."""
    data1_file_list = data1_files.split(",")
    data2_file_list = data2_files.split(",")

    # Read all benchmark files with expected count
    result1_list = [
        read_benchmark(file, n_times, n_threads, platform_override)
        if idx < len(data1_file_list)
        else None
        for idx, file in enumerate(data1_file_list)
    ]
    result2_list = [
        read_benchmark(file, n_times, n_threads, platform_override)
        if idx < len(data2_file_list)
        else None
        for idx, file in enumerate(data2_file_list)
    ]

    # Check if we have any valid data
    if not any(r and r.is_valid for r in result1_list + result2_list):
        logger.error("No valid data found in any input files")
        # Still try to generate plots with error messages

    # Generate CPU time plots
    cpu_time_fig = make_protocol_boxplots(
        result1_list,
        result2_list,
        "CPU Time Comparison",
        "CPU Time (s)",
        metric="times",
        protocols=protocols,
        use_cpu_time=True,
    )
    cpu_time_fig.savefig(f"{base_path}_cputime{extension}")
    logger.info(f"Saved CPU time plot to {base_path}_cputime{extension}")

    # Generate wall time plots
    wall_time_fig = make_protocol_boxplots(
        result1_list,
        result2_list,
        "Wall Time Comparison",
        "Wall Time (s)",
        metric="times",
        protocols=protocols,
        use_cpu_time=False,
    )
    wall_time_fig.savefig(f"{base_path}_walltime{extension}")
    logger.info(f"Saved wall time plot to {base_path}_walltime{extension}")

    # Generate memory plots (unchanged)
    mem_fig = make_protocol_boxplots(
        result1_list,
        result2_list,
        "Memory Usage Comparison",
        "Memory (MB)",
        metric="memory",
        protocols=protocols,
    )
    mem_fig.savefig(f"{base_path}_memory{extension}")
    logger.info(f"Saved memory plot to {base_path}_memory{extension}")

    # Generate efficiency metrics plots - for each protocol separately
    for idx, protocol in enumerate(protocols):
        if idx < len(result1_list) and idx < len(result2_list):
            result1 = result1_list[idx]
            result2 = result2_list[idx]

            if result1 and result1.is_valid and result2 and result2.is_valid:
                # Extract efficiency metrics
                efficiency_metrics1 = {
                    "parallelization_efficiency": result1.get_metric(
                        "parallelization_efficiency", "kernel"
                    ),
                    "memory_efficiency": result1.get_metric(
                        "memory_efficiency", "kernel"
                    ),
                    "uss_rss_ratio": result1.get_metric("uss_rss_ratio", "kernel"),
                }

                efficiency_metrics2 = {
                    "parallelization_efficiency": result2.get_metric(
                        "parallelization_efficiency", "kernel"
                    ),
                    "memory_efficiency": result2.get_metric(
                        "memory_efficiency", "kernel"
                    ),
                    "uss_rss_ratio": result2.get_metric("uss_rss_ratio", "kernel"),
                }

                # Create a protocol-specific output prefix
                protocol_base_path = f"{base_path}_{protocol}"

                # Add efficiency metrics plot
                plot_efficiency_metrics(
                    efficiency_metrics1,
                    efficiency_metrics2,
                    title=f"{protocol.upper()} Efficiency Metrics",
                    output_prefix=f"{protocol_base_path}_efficiency",
                )

                # Add workload analysis plot
                plot_workload_analysis(
                    result1,
                    result2,
                    title=f"{protocol.upper()} Workload Analysis",
                    output_prefix=protocol_base_path,
                    n_threads=n_threads,
                )


def make_protocol_boxplots(
    result1_list,
    result2_list,
    title,
    ylabel,
    metric="times",
    protocols=["chip", "rna", "wgs"],
    use_cpu_time=True,
):
    """Create boxplots for multiple protocols."""
    fig, axes = plt.subplots(1, len(protocols), figsize=(15, 5))
    fig.suptitle(title)

    # Handle case of only one protocol
    if len(protocols) == 1:
        axes = [axes]

    has_any_errors = any(r.is_error for r in result1_list + result2_list if r)
    plot_title = title
    if has_any_errors and not plot_title.endswith(" (Errors)"):
        plot_title += " (Errors)"

    fig.suptitle(plot_title)

    for idx, protocol in enumerate(protocols):
        # Get the appropriate axis
        ax = axes[idx]

        # Capitalize protocol name for display
        display_protocol = (
            protocol.upper() if len(protocol) <= 3 else protocol.capitalize()
        )

        # Handle missing data
        if idx >= len(result1_list) or idx >= len(result2_list):
            logger.warning(f"Not enough data for protocol {protocol}")
            ax.text(
                0.5,
                0.5,
                f"Missing data files for {display_protocol}",
                ha="center",
                va="center",
                fontsize=12,
                color="red",
            )
            ax.set_title(f"{display_protocol}")
            continue

        # Get the results
        result1 = result1_list[idx]
        result2 = result2_list[idx]

        # Create the boxplot with specified time type
        create_boxplot(
            ax,
            result1,
            result2,
            metric,
            title=display_protocol,
            ylabel=ylabel if idx == 0 else None,
            is_first=(idx == 0),
            use_cpu_time=use_cpu_time,
        )

    plt.tight_layout()
    return fig


def make_boxplot(
    result1,
    result2,
    title,
    ylabel,
    metric="times",
    labels=None,
    use_cpu_time=True,
):
    """
    Create a boxplot comparing two datasets.

    Args:
        result1: First BenchmarkResult
        result2: Second BenchmarkResult
        title: Plot title
        ylabel: Y-axis label
        metric: Metric to plot ("times" for time or "memory" for memory)
        labels: Labels for the boxplots (defaults to global version labels)
        use_cpu_time: If True, use CPU time; if False, use wall time
    """
    if labels is None:
        labels = [PREV_VERSION_LABEL, CURRENT_VERSION_LABEL]

    fig, ax = plt.subplots(figsize=(8, 6))

    create_boxplot(
        ax,
        result1,
        result2,
        metric,
        title,
        ylabel,
        labels=labels,
        use_cpu_time=use_cpu_time,
    )

    plt.tight_layout()
    return fig


def create_boxplot(
    ax,
    result1,
    result2,
    metric="times",
    title=None,
    ylabel=None,
    is_first=True,
    labels=None,
    use_cpu_time=True,
):
    """Helper function to create a boxplot on a given axis."""
    if labels is None:
        labels = [PREV_VERSION_LABEL, CURRENT_VERSION_LABEL]

    # Choose the appropriate metric based on what we want to show
    if metric == "times":
        if use_cpu_time:
            metric_key = "cpu_time"
            source = "kernel"
        else:
            metric_key = "wall_time"
            source = "python"  # Could also use "kernel" since both have wall_time
    elif metric == "memory":
        metric_key = "max_rss"
        source = "kernel"
    else:
        # Use the metric name directly for other metrics
        metric_key = metric
        source = "kernel"  # Default to kernel data

    data1 = result1.get_metric(metric_key, source) if result1 else []
    data2 = result2.get_metric(metric_key, source) if result2 else []

    # Flag if any result had errors
    has_errors = (result1 and result1.is_error) or (result2 and result2.is_error)

    # Check for measurement count issues
    count_issues = False
    validation_issues1 = result1.data.get("validation_issues", {}) if result1 else {}
    validation_issues2 = result2.data.get("validation_issues", {}) if result2 else {}

    if metric_key in validation_issues1 or metric_key in validation_issues2:
        count_issues = True

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
            f"No valid data{' (Error)' if result1 and result1.is_error else ''}",
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
            f"No valid data{' (Error)' if result2 and result2.is_error else ''}",
            ha="center",
            va="bottom",
            fontsize=10,
            color="red",
            transform=ax.transAxes,
        )

    # Create boxplot - use tick_labels instead of labels for matplotlib 3.9+ compatibility
    try:
        # Try the new parameter name first (matplotlib 3.9+)
        bp = ax.boxplot(
            [valid_data1, valid_data2], patch_artist=True, tick_labels=labels
        )
    except TypeError:
        # Fall back to the old parameter name
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

    # Add error/count issue indicators to title
    plot_title = title or ""
    if has_errors and count_issues:
        if not plot_title.endswith(" (Errors, Missing measurements)"):
            plot_title += " (Errors, Missing measurements)"
    elif has_errors:
        if not plot_title.endswith(" (Errors)"):
            plot_title += " (Errors)"
    elif count_issues:
        if not plot_title.endswith(" (Missing measurements)"):
            plot_title += " (Missing measurements)"

    if is_first and ylabel:
        ax.set_ylabel(ylabel)
    if plot_title:
        ax.set_title(plot_title)

    return bp


def plot_efficiency_metrics(
    metrics1, metrics2, labels=None, title=None, output_prefix=None
):
    """
    Create comparative plots for efficiency metrics

    Args:
        metrics1: Dictionary with efficiency metrics for tool 1
        metrics2: Dictionary with efficiency metrics for tool 2
        labels: Labels for the plot (defaults to global version labels)
        title: Title for the plot
        output_prefix: Prefix for output file
    """
    if labels is None:
        labels = [PREV_VERSION_LABEL, CURRENT_VERSION_LABEL]
    metrics = [
        (
            "parallelization_efficiency",
            "Parallelization Efficiency",
            "Higher values for CPU-bound tasks, lower for I/O-bound",
        ),
        (
            "memory_efficiency",
            "Memory Usage Efficiency",
            "Higher means better memory sharing across threads",
        ),
        (
            "uss_rss_ratio",
            "USS/RSS Ratio",
            "Lower means better memory sharing (private vs total memory)",
        ),
    ]

    for metric_name, metric_title, explanation in metrics:
        plt.figure(figsize=(10, 6))

        # Create boxplots
        data = [
            [v for v in metrics1[metric_name] if v is not None],
            [v for v in metrics2[metric_name] if v is not None],
        ]

        if not any(data) or all(len(d) == 0 for d in data):
            plt.text(
                0.5,
                0.5,
                "No data available for this metric",
                ha="center",
                va="center",
                fontsize=14,
                color="red",
            )
            plt.title(f"{title} - {metric_title} (No Data)")
        else:
            bp = plt.boxplot(data, labels=labels, patch_artist=True)

            # Add color coding
            for patch, color in zip(bp["boxes"], ["lightblue", "lightgreen"]):
                patch.set_facecolor(color)

            # Add individual points for transparency
            for i, d in enumerate(data):
                if d:
                    x = i + 1 + np.random.normal(0, 0.04, len(d))  # Add jitter
                    plt.scatter(x, d, alpha=0.5, color="darkblue")

            # Add mean values as text
            for i, d in enumerate(data):
                if d:
                    mean_val = np.mean(d)
                    plt.text(
                        i + 1,
                        np.max(d) * 1.05,
                        f"Mean: {mean_val:.3f}",
                        ha="center",
                        va="bottom",
                        fontweight="bold",
                    )

            # Add reference lines for parallelization efficiency
            if metric_name == "parallelization_efficiency":
                plt.axhline(
                    y=1.0,
                    color="green",
                    linestyle="--",
                    alpha=0.5,
                    label="Ideal CPU-bound",
                )
                plt.axhline(
                    y=0.5,
                    color="orange",
                    linestyle="--",
                    alpha=0.5,
                    label="Potential I/O bottleneck",
                )
                plt.legend()

            # Add explanation text
            plt.figtext(
                0.5,
                0.01,
                explanation,
                ha="center",
                fontsize=9,
                bbox={"facecolor": "lightyellow", "alpha": 0.5, "pad": 5},
            )

        plt.title(f"{title} - {metric_title}")
        plt.grid(True, linestyle="--", alpha=0.7)

        # Save plot
        output_file = f"{output_prefix}_{metric_name}.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        logger.info(f"Saved {output_file}")


def plot_workload_analysis(
    result1, result2, labels=None, title=None, output_prefix=None, n_threads=None
):
    """
    Create a plot showing workload characteristics

    Args:
        result1: First BenchmarkResult object
        result2: Second BenchmarkResult object
        labels: Labels for the two datasets (defaults to global version labels)
        title: Title for the plot
        output_prefix: Prefix for output file
        n_threads: Number of threads used (if known)
    """
    if labels is None:
        labels = [PREV_VERSION_LABEL, CURRENT_VERSION_LABEL]
    plt.figure(figsize=(10, 6))

    # Get wall and CPU times
    wall_times1 = result1.get_metric("wall_time", "python") if result1 else []
    cpu_times1 = result1.get_metric("cpu_time", "kernel") if result1 else []
    wall_times2 = result2.get_metric("wall_time", "python") if result2 else []
    cpu_times2 = result2.get_metric("cpu_time", "kernel") if result2 else []

    # Calculate ratios
    ratios1 = [
        cpu / wall
        for cpu, wall in zip(cpu_times1, wall_times1)
        if wall is not None and cpu is not None and wall > 0
    ]
    ratios2 = [
        cpu / wall
        for cpu, wall in zip(cpu_times2, wall_times2)
        if wall is not None and cpu is not None and wall > 0
    ]

    if not ratios1 and not ratios2:
        plt.text(
            0.5,
            0.5,
            "No data available for workload analysis",
            ha="center",
            va="center",
            fontsize=14,
            color="red",
        )
        plt.title(f"{title} - Workload Analysis (No Data)")
    else:
        # Plot histograms of CPU/Wall ratios
        max_ratio = max(max(ratios1, default=0), max(ratios2, default=0))
        bins = np.linspace(0, max_ratio * 1.1, 30)
        plt.hist([ratios1, ratios2], bins=bins, alpha=0.7, label=labels)

        # Add reference line for ideal CPU-bound (ratio = n_threads)
        if n_threads:
            plt.axvline(
                x=n_threads,
                color="green",
                linestyle="--",
                label=f"Ideal CPU-bound ({n_threads} threads)",
            )
            plt.axvline(
                x=n_threads * 0.7,
                color="orange",
                linestyle="--",
                label="I/O-bound threshold",
            )
            plt.axvline(
                x=n_threads * 1.3,
                color="red",
                linestyle="--",
                label="Super-linear threshold",
            )

        plt.xlabel("CPU Time / Wall Time Ratio")
        plt.ylabel("Frequency")
        plt.title(f"{title} - Workload Characteristics Analysis")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.5)

        # Add workload type annotation
        workload_types = []
        for i, ratios in enumerate([ratios1, ratios2]):
            if ratios:
                median_ratio = np.median(ratios)
                if n_threads:
                    ratio_of_expected = median_ratio / n_threads
                    if ratio_of_expected < 0.7:
                        wtype = "I/O-bound"
                    elif ratio_of_expected <= 1.3:
                        wtype = "CPU-bound"
                    else:
                        wtype = "Super-linear"

                    workload_types.append(
                        f"{labels[i]}: {wtype} (ratio: {median_ratio:.2f})"
                    )

        if workload_types:
            plt.figtext(
                0.5,
                0.01,
                " | ".join(workload_types),
                ha="center",
                fontsize=10,
                bbox={"facecolor": "lightyellow", "alpha": 0.5, "pad": 5},
            )

    # Save plot
    output_file = f"{output_prefix}_workload_analysis.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    logger.info(f"Saved {output_file}")


def get_platform(platform_override="auto"):
    """
    Get the current platform with possible override.

    Args:
        platform_override: Platform to use instead of the current one.
                          Values can be 'auto', 'Linux', or 'Darwin'

    Returns:
        String representing the platform to use for parsing
    """
    if platform_override == "auto":
        return platform.system()

    # Validate the override
    if platform_override not in ["Linux", "Darwin"]:
        logger.warning(
            f"Invalid platform override '{platform_override}', using detected platform"
        )
        return platform.system()

    logger.info(
        f"Platform override: Using '{platform_override}' instead of '{platform.system()}'"
    )
    return platform_override


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="Generate benchmark comparison plots")

    parser.add_argument(
        "--threads",
        type=int,
        required=True,
        help="Number of threads used in benchmarks",
    )

    parser.add_argument(
        "--ntimes",
        type=int,
        required=True,
        help="Number of expected measurement repetitions",
    )

    parser.add_argument(
        "--platform",
        type=str,
        choices=["auto", "Linux", "Darwin"],
        default="auto",
        help="Override platform detection for parsing results from another OS (default: auto)",
    )

    parser.add_argument(
        "output_template",
        help="Output file template (e.g., 'output/bamCoverage_human_bs100.png')",
    )
    parser.add_argument(
        "data1_files",
        help="First dataset file(s), comma-separated for multi-file comparisons",
    )
    parser.add_argument(
        "data2_files",
        help="Second dataset file(s), comma-separated for multi-file comparisons",
    )

    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_command_line_args()

    # Set logging level based on verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Get values from args
    output_template = args.output_template
    n_times = args.ntimes
    n_threads = args.threads
    platform_override = args.platform

    # Parse the output path
    output_path = Path(output_template)
    base_path = output_path.with_suffix("")

    # Log platform information
    logger.info(f"Running on: {platform.system()}")
    if platform_override != "auto":
        logger.info(f"Platform override: Using {platform_override} for parsing")

    # Process differently based on whether we're doing multi-file processing
    # (detected by comma in the input files)
    if "," in args.data1_files:
        # For bamCoverage, we know the protocols are chip,rna,wgs in that order
        protocols = ["chip", "rna", "wgs"]
        # Update these function calls to pass platform_override
        process_multiprotocol(
            args.data1_files,
            args.data2_files,
            base_path,
            protocols,
            n_times,
            n_threads,
            platform_override=platform_override,
        )
    else:
        # For single file comparisons (bamCompare, computeMatrix, multiBamSummary)
        process_single_files(
            args.data1_files,
            args.data2_files,
            base_path,
            n_times,
            n_threads,
            platform_override=platform_override,
        )

    # Rename the log file
    if cmd_being_process:
        new_log_file = log_file.replace(
            "present_results", f"present_results_{cmd_being_process}"
        )
        try:
            os.rename(log_file, new_log_file)
        except Exception as e:
            logger.warning(f"Could not rename log file: {e}")
