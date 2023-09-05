#!/bin/bash

ml load python

python join_segments_vcf.py \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta \
    ../../5-untip/unitig-mapping-3.txt \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
    utig3.fasta \
    utig1-nodes.vcf \
    utig3.vcf