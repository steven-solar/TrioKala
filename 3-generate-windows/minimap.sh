#!/bin/bash

minimap2 \
    -x asm20 \
    -o alignments/$1_$2.paf \
    alignments/fastas/$2.fasta \
    alignments/fastas/$1.fasta

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/paf_to_superset_bed.py \
    alignments/$1_$2.paf \
    alignments/$1_$2.bed
