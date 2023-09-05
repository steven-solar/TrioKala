#!/bin/bash

gunzip < Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz | awk '{if ($1==8 && $2==28730468) {print}}' > awk.txt

