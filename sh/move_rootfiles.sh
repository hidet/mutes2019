#!/bin/bash

#run=(160 161)
#for i in "${run[@]}"
for i in `seq 0 999`
do
    r=`printf %04d $i`
    #mkdir -p /Volumes/HEATES_HD/root/run$r
    mv /Volumes/SSD-PGU3/jparc2018_data/TMU_2018U/run$r/*.root /Volumes/HEATES_HD/root/run$r/
done
