#!/bin/sh

#run_n=360
run_n=412
#run_n=412
#chan_n=129
#chan_n=441
#chan_n=121

#python ana_pulse_study.py ${run_n} --chans=${chan_n} --cut_category=sec_pr_mean --cutMax=37 --energy_cut --energyMin=6900. --energyMax=6975. --pulse_ana --doiteach

#python ana_pulse_study.py ${run_n} --chans=${chan_n} --cut_category=sec_pr_mean --cutMax=37  --pulse_ana --doiteach

#python ana_pulse_study.py ${run_n} --energy_cut --energyMin=6900. --energyMax=6975 --cut_category=sec_pr_mean --cutMax=37 --doiteach --pulse_ana --chans=${chan_n}


# for group trigger max analysis by S.Y.
#python ana_pulse_study.py ${run_n} --energy_cut --energyMin=100. --energyMax=10000 --cut_category=sec_pr_mean --cutMax=37 --doiteach --pulse_ana --chans=${chan_n}
python ana_pulse_study.py ${run_n} --energy_cut --energyMin=-100. --energyMax=10000 --cut_category=sec_pr_mean --cutMax=1000 --doiteach --pulse_ana --chans=211,213,215,217

#python ana_pulse_study.py ${run_n} --energy_cut --energyMin=-100. --energyMax=10000 --cut_category=sec_pr_mean --cutMax=1000 --doiteach --pulse_ana --chans=${chan_n}
