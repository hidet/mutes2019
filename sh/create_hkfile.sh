#!/bin/sh

# 2018/06/04, S.Yamada
# cut adr_gui_log
#
# File contains
# Run_pulse,run_noise,ana_target,exttrig,grouptrig

if [ _$1 = _ ]; then
    echo "..............................................."
    echo "ERROR : need CSV file"
    echo "./create_hkfile.sh xxx.csv"
    echo "e.g. input csv file is created by "
    echo "grep -v noise csv/data_TMU_2018U.csv  > cutnoise_data_TMU_2018U.csv"
    echo "please try again"
    echo "..............................................."
    exit
fi

file=$1


ANADIR="/home/heates/ana/ana_2018_jparc"
mkdir -p ${ANADIR}/data_hk

for oneline in `cat $file`
do

echo $oneline    
   
runnum=`echo $oneline | awk -F, '{printf("%04d",$1)}'`

com="python ${ANADIR}/ana_noise.py TMU_2018U run"${runnum}
echo $com 
$com

done
