# deeptools25

> repo for benchmarking, comparison, test data and reference construction for the 25 remake.

## work-dir/

Its contents are ignored by git. Stay there.

## data/

Symlink to where all bam files, bigwig files etc are stored. Keep actual data out of git repo for size reasons. We'll make a tarball or whatever else is needed for distribution.

## scripts/

`diffbed.py` - Script to parse diff between  2 bedgraph files (new vs. old algorithm), and filter out 'false positives'. Important, when running diff, rust bedgraph should come first, and 'legacy' deeptools is second. **USAGE:** `diff rust.bedgraph original.bedgraph | python3 scripts/diffbed.py`


