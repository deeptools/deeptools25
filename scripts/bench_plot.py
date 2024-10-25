#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys

# Function to read and parse the benchmark.txt file
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

# Read data from benchmark.txt using first argument
times, memory = read_benchmark(sys.argv[1])

# Plotting
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Time plot
ax[0].bar(range(len(times)), times, color='blue')
ax[0].set_title('Average Execution Time')
ax[0].set_xlabel('Run')
ax[0].set_ylabel('Time (s)')

# Memory plot
ax[1].bar(range(len(memory)), memory, color='green')
ax[1].set_title('Memory Usage')
ax[1].set_xlabel('Run')
ax[1].set_ylabel('Memory (MB)')

plt.tight_layout()

# Save the plot using second argument
plt.savefig(sys.argv[2])
