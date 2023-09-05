#!/bin/bash

set -e

ml load python samtools minimap2 seqtk 

mkdir -p logs alignments alignments/fastas logs/alignments

compressed_paths_txt_file=compressed-nondenovo-paths.txt
compressed_paths_fasta_file=compressed-nondenovo-paths.fasta

if [[ -f $compressed_paths_txt_file && -f $compressed_paths_fasta_file ]]; then
    echo "compressed paths fasta file exists: $compressed_paths_txt_file"
else
    echo "generating compressed nondenovo paths"
    while read path_name path_nodes;
    do
        if [[ $path_name == *"_nondenovo"* ]]; then
            echo -e "$path_name\t$path_nodes"
            echo -e "$path_name\t$path_nodes" >> $compressed_paths_txt_file
        fi
    done < ../1-generate-paths/comparison-paths.txt

    python ../1-generate-paths/join_segments.py \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta \
        $compressed_paths_txt_file \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
        $compressed_paths_fasta_file
fi

compressed_nodes_txt_file=compressed-denovo-nodes.txt
compressed_nodes_fasta_file=compressed-denovo-nodes.fasta

if [[ -f $compressed_nodes_fasta_file && -f $compressed_nodes_txt_file ]]; then
    echo "compressed nodes fasta file exists: $compressed_nodes_fasta_file"
else
    echo "generating compressed denovo nodes"
    python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/get_singleton_nodes_from_paths.py \
        $compressed_paths_txt_file > $compressed_nodes_txt_file

    python ../1-generate-paths/join_segments.py \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta \
        $compressed_nodes_txt_file \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
        $compressed_nodes_fasta_file
fi

while read path_name path_nodes;
do
    fasta_file=alignments/fastas/$path_name.fasta
    if test -f $fasta_file; then
        echo "fasta file exists: $fasta_file"
    else
        touch $path_name.lst
        echo "$path_name" > $path_name.lst
        seqtk subseq \
            $compressed_paths_fasta_file \
            $path_name.lst \
            > alignments/fastas/$path_name.fasta
        rm $path_name.lst
        samtools faidx alignments/fastas/$path_name.fasta
    fi
done < $compressed_paths_txt_file

while read path_name path_nodes;
do
    fasta_file=alignments/fastas/$path_name.fasta
    if test -f $fasta_file; then
        echo "fasta file exists: $fasta_file"
    else
        touch $path_name.lst
        echo "$path_name" > $path_name.lst
        seqtk subseq \
            $compressed_nodes_fasta_file \
            $path_name.lst \
            > alignments/fastas/$path_name.fasta
        rm $path_name.lst
        samtools faidx alignments/fastas/$path_name.fasta
    fi
done < $compressed_nodes_txt_file

alignment_pairs_txt_file=alignment-pairs.txt 
python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/generate_map_pairs.py alignments/fastas > $alignment_pairs_txt_file

while read denovo non_denovo;
do
    aln_bed_file=alignments/${denovo}_$non_denovo.bed
    if test -f $aln_bed_file; then
        echo "alignment file exists: $aln_bed_file"
    else
        echo "generating alignment bed: $aln_bed_file "
        sbatch \
            --time=04:00:00 \
            --mem=121g \
            --cpus-per-task=32 \
            -o "logs/alignments/${denovo}_${non_denovo}.out" \
            minimap.sh $denovo $non_denovo
    fi
done < $alignment_pairs_txt_file

