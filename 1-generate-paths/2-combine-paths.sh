#!/bin/bash

ml load python 

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/unify_gaf_path.py \
    graph-alignment/pipeline \
    ../0-setup/multigraph.graphml \
    ../0-setup/unitig-color.csv \
    denovo-paths.txt \
    ga-paths.txt

cat bubble-paths.txt ga-paths.txt > comparison-paths.txt

cp \
    comparison-paths.txt \
    ../2-paths-consensus/6-layoutContigs/consensus_paths.txt
