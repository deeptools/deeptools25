#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

def read_benchmark(file_path):
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
            
    return {
        'times': times, 'memory': memory, 'cpu_usage': cpu_usage,
        'io_in': io_in, 'io_out': io_out, 'mean_load': mean_load,
        'cpu_time': cpu_time, 'max_uss': max_uss, 'max_pss': max_pss
    }

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

def plot_cpu_usage_time(times, cpu_usage, title="CPU Usage Over Time"):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(range(len(times)), cpu_usage, marker='o', linestyle='-', color='blue')
    ax.set_title(title)
    ax.set_xlabel("Run Number")
    ax.set_ylabel("CPU Usage (%)")
    ax.grid(True)
    return fig

def plot_io_performance(io_in, io_out, title="I/O Performance"):
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(2)
    ax.bar(x, [np.mean(io_in), np.mean(io_out)], yerr=[np.std(io_in), np.std(io_out)])
    ax.set_xticks(x)
    ax.set_xticklabels(['IO In', 'IO Out'])
    ax.set_title(title)
    ax.set_ylabel('MB')
    return fig

def plot_memory_composition(rss, uss, pss, title="Memory Usage Composition"):
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(rss))
    width = 0.25
    ax.bar(x - width, rss, width, label='RSS')
    ax.bar(x, uss, width, label='USS')
    ax.bar(x + width, pss, width, label='PSS')
    ax.set_title(title)
    ax.set_ylabel('Memory (MB)')
    ax.set_xlabel('Run Number')
    ax.legend()
    return fig

def plot_cpu_load_dist(mean_load, title="CPU Load Distribution"):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(mean_load, bins=20, color='green', alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel('Mean Load (%)')
    ax.set_ylabel('Frequency')
    return fig

def plot_efficiency(wall_time, cpu_time, title="CPU Efficiency"):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(wall_time, cpu_time, alpha=0.6)
    ax.set_title(title)
    ax.set_xlabel('Wall Time (s)')
    ax.set_ylabel('CPU Time (s)')
    # Add diagonal line for perfect efficiency
    max_val = max(max(wall_time), max(cpu_time))
    ax.plot([0, max_val], [0, max_val], 'r--', label='Perfect Scaling')
    ax.legend()
    return fig


if __name__ == '__main__':
    output_template = sys.argv[1]
    data1 = read_benchmark(sys.argv[2])
    data2 = read_benchmark(sys.argv[3])

    time_fig = make_boxplot(data1['times'], data2['times'], 'Execution Time Comparison', 'Time (s)')
    time_fig.savefig(output_template.replace('.png', '_time.png'))

    mem_fig = make_boxplot(data1['memory'], data2['memory'], 'Memory Usage Comparison', 'Memory (MB)')
    mem_fig.savefig(output_template.replace('.png', '_mem.png'))

    cpu_usage_fig = plot_cpu_usage_time(data1['times'], data1['cpu_usage'])
    cpu_usage_fig.savefig(output_template.replace('.png', '_cpu_usage.png'))

    io_fig = plot_io_performance(data1['io_in'], data1['io_out'])
    io_fig.savefig(output_template.replace('.png', '_io.png'))

    mem_comp_fig = plot_memory_composition(data1['memory'], data1['max_uss'], data1['max_pss'])
    mem_comp_fig.savefig(output_template.replace('.png', '_mem_composition.png'))

    load_fig = plot_cpu_load_dist(data1['mean_load'])
    load_fig.savefig(output_template.replace('.png', '_cpu_load.png'))

    eff_fig = plot_efficiency(data1['times'], data1['cpu_time'])
    eff_fig.savefig(output_template.replace('.png', '_efficiency.png'))
