#!/bin/bash

# TMU_2019G
export MUTESDATADIR="$MUTESHOME/data/TMU_2019G"

run=(38 40 43 44 54 61 63 66 67 76 78 84 85 87 88 97 98 100)

for r in "${run[@]}"
do
    python ipy_test.py $r -c -R
    echo $r
done  
