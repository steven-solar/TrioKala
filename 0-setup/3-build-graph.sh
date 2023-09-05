#!/bin/bash

ml load python

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/build_graph.py \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
    multigraph.graphml
