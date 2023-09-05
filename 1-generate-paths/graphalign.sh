#!/bin/bash

touch $1.lst
echo "${1}" > $1.lst
seqtk subseq \
    denovo-paths.fasta \
    $1.lst \
    > graph-alignment/pipeline/$1.fasta
rm $1.lst

echo "$1 to fasta complete"

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/get_gfa_neighborhood.py \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
    ../0-setup/multigraph.graphml \
    graph-alignment/pipeline \
    $1 5

echo "$1 neighborhood gfa complete"

GraphAligner \
    -g graph-alignment/pipeline/$1_neighborhood.gfa \
    -f graph-alignment/pipeline/$1.fasta \
    -a graph-alignment/pipeline/$1_aln.gaf \
    -x dbg \
    -t $SLURM_CPUS_PER_TASK
echo "$1 alignment complete"

# python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/aln_to_paths.py ${1}
# echo "created paths"
