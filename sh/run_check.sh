#!/bin/sh

# 2018/06/11 Hide
# changed not to use adr_gui

dir=`pwd`
ANADIR="/home/heates/ana/ana_2018_jparc"
ADRLOGDIR="${ANADIR}/log/adr"
CRYTAG="TMU_2018U"
CSVFILE="${ANADIR}/csv/data_TMU_2018U.csv"


if [ _$1 = _  ]; then
    echo "..............................................."
    echo "ERROR : need a run number"
    echo "./run_check.sh XXX"
    echo "please try again"
    echo "..............................................."
    exit
fi
run=$1


echo "run check start..."

# make symbolic links
echo "start making symbolic links..."
sh ${ANADIR}/sh/lndir_for_remotely_mounted_data.sh
ret=$?
if [ $ret -eq 1 ]; then
    echo "Error in lndir_for_remotely_mounted_data.sh"
    exit
fi
# run single ana
echo "start ana ljh files"
cd ${ANADIR}
python ${ANADIR}/ana_galen_single.py ${run}
ret=$?
if [ $ret -eq 1 ]; then
    echo "Error in ana_galen_single.py"
    exit
fi

# grptrig check
echo "start grptrig checking..."
cd ${ANADIR}
python ${ANADIR}/ana_check_grptrig.py ${CRYTAG} ${run}
ret=$?
if [ $ret -eq 1 ]; then
    echo "Error in ana_check_grptrig.py"
    exit
fi

echo "start looking at a few pixels in detail..."
cd ${ANADIR}
chan=59
python ${ANADIR}/ana_cut.py ${CRYTAG} ${run} ${chan}
chan=19
python ${ANADIR}/ana_cut.py ${CRYTAG} ${run} ${chan}
ret=$?
if [ $ret -eq 1 ]; then
    echo "Error in ana_cut.py"
    exit
fi

echo "start making temperature log..."
sh ${ANADIR}/sh/plot_adr_gui.sh ${run} ${run}
ret=$?
if [ $ret -eq 1 ]; then
    echo "Error in plot_adr_gui.sh"
    exit
fi

echo "run check end. good luck."
cd ${dir}
