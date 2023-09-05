#!/bin/bash

ml load python

python find-associations.py $1 \
    /data/Phillippy/projects/genochondromatosis/steven/data/9606.protein.aliases.v11.5.txt \
    /data/Phillippy/projects/genochondromatosis/steven/data/9606.protein.links.full.v11.5.txt
