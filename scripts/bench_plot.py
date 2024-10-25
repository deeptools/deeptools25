#!/usr/bin/env python
import matplotlib.pyplot as plt
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

def make_boxplot(data1, data2, title, ylabel):
    fig, ax = plt.subplots(figsize=(8, 5))
    bp = ax.boxplot([data1, data2], patch_artist=True, labels=['bamCoverage1', 'bamCoverage2'])
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

if __name__ == '__main__':
    output_template = sys.argv[1]
    times1, memory1 = read_benchmark(sys.argv[2])
    times2, memory2 = read_benchmark(sys.argv[3])
    time_fig = make_boxplot(times1, times2, 'Execution Time Comparison', 'Time (s)')
    time_fig.savefig(output_template.replace('.png', '_time.png'))
    mem_fig = make_boxplot(memory1, memory2, 'Memory Usage Comparison', 'Memory (MB)')
    mem_fig.savefig(output_template.replace('.png', '_mem.png'))

