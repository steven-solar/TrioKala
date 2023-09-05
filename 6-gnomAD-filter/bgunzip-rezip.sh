#!/bin/bash

ml load samtools

# bgzip -d -@$SLURM_CPUS_PER_TASK Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz

# bgzip < Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf > Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz
# tabix -p vcf Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz

# bcftools view \
#     --threads $SLURM_CPUS_PER_TASK \
#     -H \
#     Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf \
#     > Homo_sapiens-GCA_009914755.4-2022_10-gnomad.no_header.vcf

bcftools reheader -h header.txt Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf > Homo_sapiens-GCA_009914755.4-2022_10-gnomad.new_header.vcf 

bgzip \
    < Homo_sapiens-GCA_009914755.4-2022_10-gnomad.new_header.vcf \
    > Homo_sapiens-GCA_009914755.4-2022_10-gnomad.new_header.vcf.gz 

tabix -p vcf Homo_sapiens-GCA_009914755.4-2022_10-gnomad.new_header.vcf.gz