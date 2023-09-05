#!/bin/bash

ml load bcftools minimap2 python

sbatch --time=1-00:00:00 --mem=247g --cpus-per-task=32 -o "2b-trace.out" 2b-trace.sh "calls/utig1-*/utig1-*.pair.compressed.vcf" "."
sbatch --time=1-00:00:00 --mem=247g --cpus-per-task=32 -o "2b-trace-windowed.out" 2b-trace.sh "calls/utig1-*/utig1-*.pair.windowed.compressed.vcf" ".windowed."

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/combine_vcf.py \
    all-nodes.lifted.windowed.vcf \
    "calls/utig1-*/utig1-*.pair.windowed.lifted.vcf"
bcftools sort all-nodes.lifted.windowed.vcf -o all-nodes.lifted.windowed.sorted.vcf
rm all-nodes.lifted.windowed.vcf

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/combine_vcf.py \
    all-nodes.lifted.vcf \
    "calls/utig1-*/utig1-*.pair.lifted.vcf"
bcftools sort all-nodes.lifted.vcf -o all-nodes.lifted.sorted.vcf
rm all-nodes.lifted.vcf

python /data/solarsj/genochondromatosis/steven/python_scripts/combine_vcf.py \
    all-nodes.lifted.SVs.vcf \
    "calls/utig1-*/utig1-*.SVs.lifted.vcf"
bcftools sort all-nodes.lifted.SVs.vcf -o all-nodes.lifted.SVs.sorted.vcf
rm all-nodes.lifted.SVs.vcf
