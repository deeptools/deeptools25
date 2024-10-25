# deeptools25

> repo for benchmarking, comparison, test data and reference construction for the 25 remake.

## scripts/

`diffbed.py` - Script to parse diff between  2 bedgraph files (new vs. old algorithm), and filter out 'false positives'. Important, when running diff, rust bedgraph should come first, and 'legacy' deeptools is second. **USAGE:** `diff rust.bedgraph original.bedgraph | python3 scripts/diffbed.py`


