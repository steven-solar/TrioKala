#!/bin/bash

set -e 
ml load python seqtk
mkdir -p lifts logs/lifts lifts/fastas lifts/chains

alignment_pairs_txt_file=alignment-pairs.txt

while read denovo_path_name nondenovo_path_name;
do
    echo $nondenovo_path_name
    touch $nondenovo_path_name.lst
    echo "$nondenovo_path_name" > $nondenovo_path_name.lst
    seqtk subseq \
        ../2-paths-consensus/7-unitigConsensus/assembly.fasta \
        $nondenovo_path_name.lst \
        > lifts/fastas/$nondenovo_path_name.original.fasta
    rm $nondenovo_path_name.lst 

    bash /data/solarsj/hp_compression_mapping/build-chains.sh \
        lifts/fastas/$nondenovo_path_name.original.fasta \
        lifts/fastas/$nondenovo_path_name.compressed.fasta \
        lifts/chains/$nondenovo_path_name
    # python /data/solarsj/hp_compression_mapping/build_sparse_compression_map.py \
    #     -c \
    #     -d 1 \
    #     lifts/fastas/$nondenovo_path_name.original.fasta \
    #     lifts/fastas/$nondenovo_path_name.compressed.fasta \
    #     lifts/chains/$nondenovo_path_name.json

    # for file in alignments/*_$nondenovo_path_name*.bed;
    # do
    #     echo $file
    #     bash /data/solarsj/hp_compression_mapping/lift.sh \
    #         $file \
    #         lifts/chains/${nondenovo_path_name}_uncompression.chn \
    #         lifts/fastas/$nondenovo_path_name.original.fasta \
    #         lifts/chains/$nondenovo_path_name \
    #         cu
        # python /data/solarsj/hp_compression_mapping/lift_seqs.py \
        #     -cu \
        #     lifts/fastas/$nondenovo_path_name.original.fasta \
        #     lifts/chains/$nondenovo_path_name.json \
        #     $file \
        #     $file.lifted
    # done
    bash /data/solarsj/hp_compression_mapping/lift.sh \
        alignments/${denovo_path_name}_$nondenovo_path_name.bed \
        lifts/chains/${nondenovo_path_name}_uncompression.chn \
        lifts/fastas/$nondenovo_path_name.original.fasta \
        lifts/chains/$nondenovo_path_name \
        cu
done < $alignment_pairs_txt_file  

