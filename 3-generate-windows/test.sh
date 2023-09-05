#!/bin/bash

ml load python ucsc

nondenovo_path_name=utig1-138458_nondenovo
for file in alignments/*_$nondenovo_path_name*.bed;
    do
        echo $file
        bash /data/solarsj/hp_compression_mapping/lift.sh \
            $file \
            lifts/chains/${nondenovo_path_name}_uncompression.chn \
            lifts/fastas/$nondenovo_path_name.original.fasta \
            lifts/chains/$nondenovo_path_name \
            cu
    done