#!/bin/bash

ml load samtools

bcftools view -i "MIN(FMT/DP)>$1" vep.vcf.gnomad_keep > final.vcf
