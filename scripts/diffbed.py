#!/usr/bin/env python3
import sys
import re


def process_diff(diff_input):
    current_diff = []
    output_diff = []
    for line in diff_input:
        if re.match(r"^\d", line):
            # If the line starts with a number, it's a new diff unit
            if current_diff:
                if not is_false_positive(current_diff):
                    output_diff.extend(current_diff)
                current_diff = []
        current_diff.append(line)

    # Check the last diff unit
    if current_diff and not is_false_positive(current_diff):
        output_diff.extend(current_diff)

    return output_diff


def is_false_positive(diff_unit):
    original_line = None
    modified_lines = []

    for line in diff_unit:
        if line.startswith("<"):
            original_line = line[2:].strip()
        elif line.startswith(">"):
            modified_lines.append(line[2:].strip())
        elif line.startswith("---"):
            continue

    if not original_line or not modified_lines:
        return False

    # Parse the original line
    orig_chrom, orig_start, orig_end, orig_value = original_line.split()
    orig_start, orig_end, orig_value = int(orig_start), int(orig_end), int(orig_value)

    # Check if all modified lines have the same chromosome and value
    first_mod = modified_lines[0].split()
    first_chrom = first_mod[0]
    first_value = int(first_mod[3])  # Convert to int

    # Check all modified lines match chromosome and value
    if any(
        line.split()[0] != first_chrom
        or int(line.split()[3]) != first_value  # Convert to int
        for line in modified_lines
    ):
        return False

    if first_chrom != orig_chrom or first_value != orig_value:
        return False

    # Check if modified lines form a continuous interval matching the original
    mod_intervals = [
        (int(line.split()[1]), int(line.split()[2])) for line in modified_lines
    ]
    mod_intervals.sort()

    # Check if intervals are continuous and match original bounds
    if mod_intervals[0][0] != orig_start:
        return False

    for i in range(len(mod_intervals) - 1):
        if mod_intervals[i][1] != mod_intervals[i + 1][0]:
            return False

    if mod_intervals[-1][1] != orig_end:
        return False

    return True


def main():
    diff_input = sys.stdin.readlines()
    filtered_diff = process_diff(diff_input)
    for line in filtered_diff:
        print(line, end="")


if __name__ == "__main__":
    main()
