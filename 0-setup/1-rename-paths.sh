#!/bin/bash

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/replacepaths_nodeID_utigID.py \
    ../../2-processGraph/unitig-mapping-1.txt \
    ../../1-buildGraph/paths.gaf \
    unitig-paths.gaf
