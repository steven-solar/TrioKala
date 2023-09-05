#!/bin/bash

ml load python

# while read sample;
# do
#     echo $sample
#     python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/get_readnames.py $sample
# done < "samples.txt"

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/get_readnames.py $1