#!/bin/bash

ml load python samtools

# bgzip \
#     < /data/Phillippy/projects/genochondromatosis/steven/verkko/assembly_trio/analysis/compare_haps/vep-2/vep.shortchrs.vcf \
#     > /data/Phillippy/projects/genochondromatosis/steven/verkko/assembly_trio/analysis/compare_haps/vep-2/vep.shortchrs.vcf.bgz

# tabix -p \
#     vcf /data/Phillippy/projects/genochondromatosis/steven/verkko/assembly_trio/analysis/compare_haps/vep-2/vep.shortchrs.vcf.bgz

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/merge_vep_gnomad.py \
    /data/Phillippy/projects/genochondromatosis/steven/data/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.new_header.vcf.gz \
    ../5-vep/all-nodes.vep.vcf 0
mv ../5-vep/all-nodes.vep.vcf.gnomad_* .
# bcftools view -i "MIN(FMT/DP)>$1" all-nodes.vep.vcf.gnomad_keep > all-nodes.final.vcf

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/merge_vep_gnomad.py \
    /data/Phillippy/projects/genochondromatosis/steven/data/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.new_header.vcf.gz \
    ../5-vep/all-nodes.windowed.vep.vcf 0
mv ../5-vep/all-nodes.windowed.vep.vcf.gnomad_* .
# bcftools view -i "MIN(FMT/DP)>$1" all-nodes.vep.vcf.gnomad_keep > all-nodes.final.vcf
