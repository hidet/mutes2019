#!/bin/sh

# 2018/06/08, H.Noda
# copy adr_gui_log from pcuser to nuc1
# 2018/06/16, Hide modified ....

LOGDIR="/home/heates/ana/ana_2018_jparc/log/adr"
rsync -arvu pcuser@10.105.51.249:/home/pcuser/adr_gui/adr_gui_log ${LOGDIR}/

day0=`date -d yesterday | awk '{if ($3/10 < 1) print "0"$3; else print $3}'`
day1=`date | awk '{if ($3/10 < 1) print "0"$3; else print $3}'`

grep 201806${day0} ${LOGDIR}/adr_gui_log > ${LOGDIR}/201806${day0}_adr_gui_log
grep 201806${day1} ${LOGDIR}/adr_gui_log > ${LOGDIR}/201806${day1}_adr_gui_log
cat ${LOGDIR}/201806${day0}_adr_gui_log ${LOGDIR}/201806${day1}_adr_gui_log > ${LOGDIR}/201806${day0}-${day1}_adr_gui_log
