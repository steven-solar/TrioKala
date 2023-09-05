#!/bin/bash

set -e

ml load python graphaligner seqtk

mkdir -p graph-alignment graph-alignment/pipeline

graph_fasta_file=../../2-processGraph/unitig-unrolled-hifi-resolved.fasta

if test -f $graph_fasta_file; then 
    echo "graph fasta file exists: $graph_fasta_file"
else
    python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/gfa_to_fasta.py \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta
fi 

denovo_paths_file=denovo-paths.txt
if test -f $denovo_paths_file; then
    echo "denovo paths file exists: $denovo_paths_file"
else
    python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/expand_denovo_paths.py \
        ../0-setup/unitig-color.csv \
        ../0-setup/unitig-denovos.filtered.txt \
        ../0-setup/multigraph.graphml \
        denovo-paths.txt
fi

denovo_fasta_file=denovo-paths.fasta
if test -f $denovo_fasta_file; then
    echo "denovo fasta file exists: $denovo_fasta_file"
else
    python join_segments.py \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta \
        denovo-paths.txt \
        ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
        denovo-paths.fasta
fi

bubble_paths_file=bubble-paths.txt
if test -f $bubble_paths_file; then
    echo "bubble paths file exists: $bubble_paths_file"
else
    python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/find_simple_bubbles.py \
        ../0-setup/unitig-color.csv \
        denovo-paths.txt \
        ../0-setup/multigraph.graphml \
        bubble-paths.txt
fi

while read DENOVO_NODE DENOVO_PATH
do
    if grep -q $DENOVO_PATH bubble-paths.txt; then
        echo "$DENOVO_NODE    $DENOVO_PATH already found in bubble-paths.txt"
    else
        aln_file=graph-alignment/pipeline/${DENOVO_NODE}_aln.gaf
        if test -f $aln_file; then
            echo "$DENOVO_NODE already GraphAligned: $aln_file"
        else
            echo "$DENOVO_NODE being GraphAligned"
            sbatch \
                --time=04:00:00 \
                --partition=norm,quick \
                --mem=121g \
                --cpus-per-task=32 \
                -o "graph-alignment/pipeline/$DENOVO_NODE.out" \
                graphalign.sh $DENOVO_NODE
        fi
    fi
done < denovo-paths.txt
