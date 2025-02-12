#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import platform
import sys

def read_benchmark(file_path):
    command = file_path.split('/')[1].split('_')[0][:-1]  # this could've been a RegEx :/
    times = []
    memory = []
    cpu_usage = []
    io_in = []
    io_out = []
    mean_load = []
    cpu_time = []
    max_uss = []
    max_pss = []
    
    with open(file_path, 'r') as file:
        header = next(file)
        for line in file:
            parts = line.strip().split('\t')
            times.append(float(parts[0]))
            memory.append(float(parts[2]))
            max_uss.append(float(parts[4]))
            max_pss.append(float(parts[5]))
            io_in.append(float(parts[6]))
            io_out.append(float(parts[7]))
            mean_load.append(float(parts[8]))
            cpu_time.append(float(parts[9]))
            cpu_usage.append(float(parts[15].split()[0]))

    # This is a workaround for macOS, where the memory values are all zeroes...
    # although Snakemake uses psutil, and it claims to support macOS :/
    # https://psutil.readthedocs.io/en/latest/index.html#psutil.Process.memory_info
    if platform.system() == "Darwin":
        Ntimes = len(memory)
        assert sum(memory) == 0, "Memory values are not all zeroes"
        log_files = glob.glob(f"logs/{command}[1-2]_[0-9]+.txt")
        memory = parse_memory_from_logs(log_files)
        assert len(memory) == Ntimes, "Expected {} memory values, got {}".format(Ntimes, len(memory))

    return {
        'times': times, 'memory': memory, 'cpu_usage': cpu_usage,
        'io_in': io_in, 'io_out': io_out, 'mean_load': mean_load,
        'cpu_time': cpu_time, 'max_uss': max_uss, 'max_pss': max_pss
    }

def parse_memory_from_logs(log_files):
    memory_values = []
    for log_file in log_files:
        with open(log_file, 'r') as file:
            for line in file:
                if "maximum resident set size" in line:
                    memory_value = int(line.split()[0])
                    memory_values.append(memory_value)
                    break  # Go to next log file
    return memory_values

def create_benchmark_df(file_paths_python, file_paths_rust):
    dfs = []
    for cmd in ['bamCoverage', 'bamCompare', 'computeMatrix', 'multiBamSummary']:
        data_py = read_benchmark(file_paths_python[cmd])
        df_py = pd.DataFrame({
            'Cmd': cmd,
            'Backend': 'Python',
            'Trial': range(1, len(data_py['times']) + 1),
            'Time': data_py['times'],
            'Memory': data_py['memory']
        })
        data_rust = read_benchmark(file_paths_rust[cmd])
        df_rust = pd.DataFrame({
            'Cmd': cmd,
            'Backend': 'Rust',
            'Trial': range(1, len(data_rust['times']) + 1),
            'Time': data_rust['times'],
            'Memory': data_rust['memory']
        })
        dfs.extend([df_py, df_rust])
    return pd.concat(dfs, ignore_index=True)

def make_boxplot(data1, data2, title, ylabel):
    fig, ax = plt.subplots(figsize=(8, 5))
    bp = ax.boxplot([data1, data2], patch_artist=True, labels=['Legacy', 'Maturin'])
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    colors = ['lightblue', 'lightgreen']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    plt.setp(bp['medians'], color='navy')
    plt.setp(bp['whiskers'], color='navy')
    plt.setp(bp['caps'], color='navy')
    plt.setp(bp['fliers'], marker='o', markerfacecolor='red', alpha=0.5)
    return fig

def make_protocol_boxplots(data1_files, data2_files, title, ylabel, metric='times', protocols=['ChIP', 'RNA', 'WGS']):
    fig, axes = plt.subplots(1, len(protocols), figsize=(15, 5))
    fig.suptitle(title)
    
    for idx, (protocol, d1_file, d2_file) in enumerate(zip(protocols, data1_files.split(','), data2_files.split(','))):
        data1 = read_benchmark(d1_file)
        data2 = read_benchmark(d2_file)
        
        bp = axes[idx].boxplot([data1[metric], data2[metric]], patch_artist=True, labels=['Legacy', 'Maturin'])
        axes[idx].set_title(f'{protocol}')
        axes[idx].set_ylabel(ylabel if idx == 0 else '')
        
        colors = ['lightblue', 'lightgreen']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        plt.setp(bp['medians'], color='navy')
        plt.setp(bp['whiskers'], color='navy')
        plt.setp(bp['caps'], color='navy')
        plt.setp(bp['fliers'], marker='o', markerfacecolor='red', alpha=0.5)
    
    plt.tight_layout()
    return fig


if __name__ == '__main__':
    output_template = sys.argv[1]
    if ',' in sys.argv[2]:
        time_fig = make_protocol_boxplots(sys.argv[2], sys.argv[3], 'Execution Time Comparison', 'Time (s)', metric='times')
        time_fig.savefig(output_template.replace('.png', '_time.png'))
        mem_fig = make_protocol_boxplots(sys.argv[2], sys.argv[3], 'Memory Usage Comparison', 'Memory (MB)', metric='memory')
        mem_fig.savefig(output_template.replace('.png', '_mem.png'))
    else:
        data1 = read_benchmark(sys.argv[2])
        data2 = read_benchmark(sys.argv[3])
        
        time_fig = make_boxplot(data1['times'], data2['times'], 'Execution Time Comparison', 'Time (s)')
        time_fig.savefig(output_template.replace('.png', '_time.png'))
        
        mem_fig = make_boxplot(data1['memory'], data2['memory'], 'Memory Usage Comparison', 'Memory (MB)')
        mem_fig.savefig(output_template.replace('.png', '_mem.png'))
    
    # dfs = []
    # for cmd_pair in zip(['bamCoverage', 'bamCompare', 'computeMatrix', 'multiBamSummary'], 
    #                     sys.argv[2].split(','), 
    #                     sys.argv[3].split(',')):
    #     cmd, py_file, rust_file = cmd_pair
    #     dfs.append(create_benchmark_df(py_file, cmd, 'Python', 1))
    #     dfs.append(create_benchmark_df(py_file, cmd, 'Python', 2))
    #     dfs.append(create_benchmark_df(rust_file, cmd, 'Rust', 1))
    #     dfs.append(create_benchmark_df(rust_file, cmd, 'Rust', 2))

    # final_df = pd.concat(dfs, ignore_index=True)
    # final_df.to_csv('benchmark_results.csv', index=False)
