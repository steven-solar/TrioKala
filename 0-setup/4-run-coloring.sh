#!/usr/bin/sh

readnames=""
while read sample;
do
    readnames+="$sample.readnames.txt "
done < samples.txt

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/color_graph.py \
    $readnames \
    < unitig-paths.gaf \
    > unitig-color.csv
