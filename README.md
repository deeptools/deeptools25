# deeptools25

> repo for benchmarking, comparison, test data and reference construction for the 25 remake.

## Install

1. _(mac only: `xcode-select --install`)_
1. `sudo apt install clang`
1. create and activate a conda env with `python==3.9.12 deeptools==3.5.5 pybigwig==0.3.22`
1. `git clone git@github.com:WardDeb/deepTools.git deepToolsWard && cd deepToolsWard`
1. `git checkout -b maturin && pip install -e .`

## Contents

#### work-dir/

Its contents are ignored by git. Stay there.

#### data/

Symlink to where all bam files, bigwig files etc are stored. Keep actual data out of git repo for size reasons. We'll make a tarball or whatever else is needed for distribution.

#### scripts/

`diffbed.py` - Script to parse diff between  2 bedgraph files (new vs. old algorithm), and filter out 'false positives'. Important, when running diff, rust bedgraph should come first, and 'legacy' deeptools is second. **USAGE:** `diff rust.bedgraph original.bedgraph | python3 scripts/diffbed.py`


