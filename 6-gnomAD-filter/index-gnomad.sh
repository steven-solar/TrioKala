#!/bin/bash

ml load samtools 

cd /data/Phillippy/projects/genochondromatosis/steven/data
tabix -p vcf Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz

