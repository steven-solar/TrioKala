#!/bin/bash 

# set -e 

ml load samtools minimap2 seqtk cmake ucsc python
source myconda
mamba activate ccs_project

mkdir -p fastas logs calls

# java -jar ~/picard.jar CreateSequenceDictionary R=/data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta O=/data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.dict

while read denovo_name nondenovo_name;
do
    sbatch \
        --time=1-00:00:00 \
        --mem=121g \
        --cpus-per-task=32 \
        -o "logs/${denovo_name}_$nondenovo_name.out" \
        1b-run.sh $denovo_name $nondenovo_name
    # bash 1b-run.sh $denovo_name $nondenovo_name 2
done < ../3-generate-windows/alignment-pairs.txt
