#!/bin/sh

ANADIR="/home/heates/ana/ana_2018_jparc"
LOGDIR="${ANADIR}/log/adr"
CRYTAG="TMU_2018U"
CSVFILE="${ANADIR}/csv/data_TMU_2018U.csv"

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

cd ${ANADIR}
# -----------------------------------
rsync -arvu pcuser@10.105.51.249:/home/pcuser/adr_gui/adr_gui_log ${LOGDIR}/
day0=`date -d yesterday | awk '{if ($3/10 < 1) print "0"$3; else print $3}'`
day1=`date | awk '{if ($3/10 < 1) print "0"$3; else print $3}'`
grep 201806${day0} ${LOGDIR}/adr_gui_log > ${LOGDIR}/201806${day0}_adr_gui_log
grep 201806${day1} ${LOGDIR}/adr_gui_log > ${LOGDIR}/201806${day1}_adr_gui_log
cat ${LOGDIR}/201806${day0}_adr_gui_log ${LOGDIR}/201806${day1}_adr_gui_log > ${LOGDIR}/201806${day0}-${day1}_adr_gui_log
# -----------------------------------
grep -v noise ${CSVFILE} | awk -F"," '($1 > '${runstart}'-1 && $1 < '${runend}'+1){print $1","$2","$3","$4","$5}' > ${ANADIR}/data_hk/cutnoise_cutold_data_TMU_2018U.csv
sh ${ANADIR}/sh/create_hkfile.sh ${ANADIR}/data_hk/cutnoise_cutold_data_TMU_2018U.csv
cd ${ANADIR}/data_hk
awk '(NR%10==0){print $0}' *ch135* > cat_per10_ch135.csv
awk '(NR%10==0){print $0}' *ch225* > cat_per10_ch225.csv
awk '(NR%10==0){print $0}' *ch323* > cat_per10_ch323.csv
awk '(NR%10==0){print $0}' *ch43* > cat_per10_ch43.csv
ls ${ANADIR}/data_hk/cat_per10_ch* > ${ANADIR}/data_hk/p_rms_list_per10.list
cd ${ANADIR}
python ${ANADIR}/plot_adrgui.py ${LOGDIR}/201806${day0}-${day1}_adr_gui_log ${CSVFILE} ${ANADIR}/data_hk/p_rms_list_per10.list
# -----------------------------------
