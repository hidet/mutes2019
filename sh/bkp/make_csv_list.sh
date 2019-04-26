#!/bin/sh

# 2018/06/08, H.Noda
# Take csv file, remove noise run, create HK file, and make pretrigger list

# 2018/06/11 Hide fixed the bug ($2 is noise run)
# 2018/06/16 Hide brushup

dir=`pwd`
ANADIR="/home/heates/ana/ana_2018_jparc"

if [ _$1 = _  ]; then
    echo "..............................................."
    echo "ERROR : need the start run number today"
    echo "./create_hkfile.sh (the start run number today)"
    echo "please try again"
    echo "..............................................."
    exit

elif [ _$2 = _ ]; then
    echo "..............................................."
    echo "ERROR : need the end run number today"
    echo "./create_hkfile.sh (the end run number today)"
    echo "please try again"
    echo "..............................................."
    exit
fi

runstart=$1
runend=$2
csvfile="/home/heates/ana/ana_2018_jparc/csv/data_TMU_2018U.csv"
cd ${ANADIR}
#grep -v noise ${csvfile} | awk -F"," '($1 > '${runstart}'-1 && $1 < '${runend}'+1){print $1","$2","$3","$4","$5}' > ${ANADIR}/cutnoise_cutold_data_TMU_2018U.csv
#./create_hkfile.sh cutnoise_cutold_data_TMU_2018U.csv
grep -v noise ${csvfile} | awk -F"," '($1 > '${runstart}'-1 && $1 < '${runend}'+1){print $1","$2","$3","$4","$5}' > ${ANADIR}/data_hk/cutnoise_cutold_data_TMU_2018U.csv
sh ${ANADIR}/sh/create_hkfile.sh ${ANADIR}/data_hk/cutnoise_cutold_data_TMU_2018U.csv

cd ${ANADIR}/data_hk
awk '(NR%10==0){print $0}' *ch135* > cat_per10_ch135.csv
awk '(NR%10==0){print $0}' *ch225* > cat_per10_ch225.csv
awk '(NR%10==0){print $0}' *ch323* > cat_per10_ch323.csv
awk '(NR%10==0){print $0}' *ch43* > cat_per10_ch43.csv
ls ${ANADIR}/data_hk/cat_per10_ch* > ${ANADIR}/data_hk/p_rms_list_per10.list
cd ${dir}
