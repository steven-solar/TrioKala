#!/bin/bash

ml load python

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/filter_singletons_tips.py \
    unitig-color.csv \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
    denovos.filtered.txt \
    denovos.tips.txt \
    denovos.singletons.txt
