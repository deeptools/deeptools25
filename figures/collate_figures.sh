#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

## Figure 1 - overview
cd "$SCRIPT_DIR/figure1"
python fig1_subfigs.py
typst compile -f png --ppi 300 fig1.typ "$SCRIPT_DIR/figure1.png"
typst compile -f pdf --ppi 300 fig1.typ "$SCRIPT_DIR/figure1.pdf"

## Figure 2 - benchmark (generated before, results committed to repo)
cp "$SCRIPT_DIR/../benchmark/results/performance.png" "$SCRIPT_DIR/figure2.png"
cp "$SCRIPT_DIR/../benchmark/results/performance.pdf" "$SCRIPT_DIR/figure2.pdf"

## Figure 3 - galaxy, to implement

## Figure 4 - example (generated before, results committed to repo)cp
cp "$SCRIPT_DIR/../example/results/figure_combined.png" "$SCRIPT_DIR/figure4.png"
cp "$SCRIPT_DIR/../example/results/figure_combined.pdf" "$SCRIPT_DIR/figure4.pdf"
