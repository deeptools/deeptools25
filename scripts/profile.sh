#!/usr/bin/env bash
PATHROOT=/data/manke/processing/salatino/deepdev/deeptools25/work-dir/testfiles
samply record bamCoverage2 -b ${PATHROOT}/S2_H2AZ.bam \
    -o ${PATHROOT}/output/new.bg -of bedgraph -bs 1 -p 1
