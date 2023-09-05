#!/bin/bash

ml load python

python path_to_node_vcf.py \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta \
    ../1-generate-paths/comparison-paths.txt \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
    path-to-node.fasta \
    all.pair.windowed.compressed.vcf \
    utig1-nodes.vcf