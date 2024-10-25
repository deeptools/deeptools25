#!/usr/bin/env python
import matplotlib.pyplot as plt
# import numpy as np
import sys


def read_benchmark(file_path):
    times = []
    memory = []
    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            parts = line.split()
            times.append(float(parts[0]))  # 's' column for times
            memory.append(float(parts[2]))  # 'max_rss' column for memory usage
    return times, memory


def make_boxplot(data, title, ylabel):
    fig, ax = plt.subplots(figsize=(6, 5))
    bp = ax.boxplot([data], patch_artist=True)
    # ax.plot([1], [np.mean(data)], marker='*', color='red', markersize=10, label='Mean')
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    # ax.legend()
    plt.setp(bp['boxes'], facecolor='lightblue', alpha=0.7)
    plt.setp(bp['medians'], color='navy')
    plt.setp(bp['whiskers'], color='navy')
    plt.setp(bp['caps'], color='navy')
    plt.setp(bp['fliers'], marker='o', markerfacecolor='red', alpha=0.5)
    return fig

if __name__ == '__main__':
    times, memory = read_benchmark(sys.argv[1])
    output_template = sys.argv[2]

    time_fig = make_boxplot(times, 'Execution Time Distribution', 'Time (s)')
    time_fig.savefig(output_template.replace('.png', '_time.png'))

    mem_fig = make_boxplot(memory, 'Memory Usage Distribution', 'Memory (MB)')
    mem_fig.savefig(output_template.replace('.png', '_mem.png'))

