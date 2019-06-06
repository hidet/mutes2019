#!/bin/bash

# TMU_2019G
export MUTESDATADIR="$MUTESHOME/data/TMU_2019G"

#python ipy_test.py 36 -ceg -REG --beam=off --sprmc=on --jbrsc=on
#python ipy_test.py 58,59 -ceg -REG --beam=off --sprmc=on --jbrsc=on
#python ipy_test.py 79,80,81 -ceg -REG --beam=off --sprmc=on --jbrsc=on
python ipy_test.py 79,80,81 -fscegd -REG --beam=off --sprmc=on --jbrsc=on
