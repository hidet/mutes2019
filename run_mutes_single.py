#!/usr/bin/env python

""" run_mutes_single.py is a basic analysis tools for MUTES project. 

History: 
2018-07-04 ; ver 1.0; made from ana_galen_single.py, S.Yamada
2018-07-05 ; ver 1.1; update to dump ROOT and pulses 
2018-07-06 ; ver 1.2; a minor update
2019-01-31 ; ver 1.3; H.Tatsuno modified drastically, removed all plots functions and simplified
2019-04-26 ; ver 1.4; H.Tatsuno added multiple runs analysis
2019-05-09 ; ver 1.5; H.Tatsuno minor change
2019-05-22 ; ver 1.6; H.Tatsuno minor change

"""

__author__ =  'a bad boy HT'
__version__ = '1.6'

import matplotlib
#matplotlib.use("Agg")
import mass
import monkeypatch # patches mass to avoid crashes
import numpy as np
import pylab as plt
params = {'xtick.labelsize': 10, # x ticks
          'ytick.labelsize': 10, # y ticks
          'legend.fontsize': 8
                    }
plt.rcParams['font.family'] = 'serif' # 
plt.rcParams.update(params)
import os
import sys
import h5py
import datetime
import pandas as pd
import optparse

import mutes_ana as MUTES
import tesmap_forgrptrig as tesmap
tesmap = reload(tesmap)
MUTES = reload(MUTES) # for ipython to reload when it changed


BADCHS = [3,9,39,77,83,85,111,337,367,375,423]# initially disconnected
BADCHS.extend([117,203,233])# bad channels
BADCHS.extend([5,177,257,265,293])# strange channels
BADCHS.sort()

maxchans = 240
#maxchans = 10

usage = u'(1) %prog 278 (basic analysis), (2) %prog 278'
version = __version__
parser = optparse.OptionParser(usage=usage, version=version)
# setting for options
# how to use
# -eg -REG : with exttrig, grouptrig, and ROOT, ROOTEXT, ROOTGRP (without forceNew, summaryNew, calibNew)
# -fsceg -REG : all analysis with forceNew, summaryNew, calibNew, exttrig, grouptrig, and ROOT, ROOTEXT, ROOTGRP
parser.add_option('-f', '--force',    dest='forceNew',   action='store_true',  help='True to update filter (default=False)',      metavar='FORCE',   default=False)
parser.add_option('-s', '--summary',  dest='summaryNew', action='store_true',  help='True to update summary (default=False)',     metavar='SUMMARY', default=False)
parser.add_option('-c', '--calib',    dest='calibNew',   action='store_true',  help='True to update calibration (default=False)', metavar='CALIB',   default=False)
parser.add_option('-e', '--exttrig',  dest='externTrig', action='store_true',  help='True for calc externTrig (default=False)',   metavar='EXTTRIG', default=False)
parser.add_option('-g', '--grptrig',  dest='groupTrig',  action='store_true',  help='True for calc groupTrig (default=False)',    metavar='GRPTRIG', default=False)
parser.add_option('-R', '--dumproot', dest='dumproot',   action="store_true",  help='dump ROOT except for pulses (default=False)',metavar='DUMPROOT',default=False)
parser.add_option('-E', '--rootext',  dest='rootext',    action='store_true',  help='True ROOT with externTrig (default=False)',  metavar='ROOTEXT', default=False)
parser.add_option('-G', '--rootgrp',  dest='rootgrp',    action='store_true',  help='True ROOT with groupTrig (default=False)',   metavar='ROOTGRP', default=False)
parser.add_option('-d', '--delete',   dest='delete',     action="store_true",  help='True to delete hdf5 file (default=False)',   metavar='DELETE',  default=False)
parser.add_option('--beam',   dest='beam',     action="store",type=str, help='set beam catecut (default=None, on or off)',default="None")
parser.add_option('--sprmc',  dest='sprmc',    action="store",type=str, help='set sprmc catecut (default=None, on or off)',default="None")
parser.add_option('--jbrsc',  dest='jbrsc',    action="store",type=str, help='set jbrsc catecut (default=None, on or off)',default="None")
parser.add_option('--pre',    dest='cut_pre',  action="store",type=int, help='set cut for pre samples',default=0)
parser.add_option('--post',   dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)
parser.add_option('--adr',    dest='adr',      action="store",type=str, help='set adr tag (default=TMU_2019, TMU_2018,TMU_2019,...)',default="TMU_2019")
parser.add_option('--cool',   dest='cool',     action="store",type=str, help='set cooling tag (default=G, G,H,I,...)',default="G")

options,args = parser.parse_args()

#### get options ####
FORCE          = options.forceNew
SUMMARY        = options.summaryNew
CALIB          = options.calibNew
EXTTRIG        = options.externTrig
GRPTRIG        = options.groupTrig
DELETE         = options.delete
DUMPROOT       = options.dumproot
ROOTEXT        = options.rootext
ROOTGRP        = options.rootgrp
beam           = options.beam
sprmc          = options.sprmc
jbrsc          = options.jbrsc
cut_pre        = options.cut_pre
cut_post       = options.cut_post
adr            = options.adr
cool           = options.cool

print "[START] " + __file__
ANADIR=os.environ.get("MUTESANADIR","")
DATADIR=os.environ.get("MUTESDATADIR","")
if ANADIR == "" or DATADIR == "":
    print "[ERROR] Set MUTESANADIR and MUTESDATADIR"
    sys.exit()
tmpdd=DATADIR.split("/TMU_")
if len(tmpdd)>1: DATADIR=tmpdd[0]# if you set the enf of DATADIR as /TMU_2019X
RUNTAG=adr+cool# TMU_2019X
DATADIR=DATADIR+"/"+RUNTAG
print "RUNTAG  = ", RUNTAG
print "ANADIR  = ", ANADIR
print "DATADIR = ", DATADIR

if os.path.isdir(DATADIR)==False: 
    print "%s is missing"%DATADIR
    print "please check the options --adr=TMU_2019 or --cool=G,H,I,... "
    sys.exit()

RUNINFO="./csv/data_%s.csv"%(RUNTAG)
if os.path.exists(RUNINFO)==False: 
    print "%s is missing"%RUNINFO
    sys.exit(0)

#GRTINFO="./csv/grptrig_twocol.txt"
GRTINFO="./csv/grptrig_wo_neighbor.csv"
if os.path.exists(GRTINFO)==False: 
    print "%s is missing"%GRTINFO
    sys.exit(0)
    
COLUMN_INFO = "./csv/column_info.csv"
if os.path.exists(COLUMN_INFO)==False:
    print "%s is missing"%COLUMN_INFO
    sys.exit()


catecut = {}
prime = "on"; catecut["prime"] = prime;
if not sprmc=="None": catecut["sprmc"] = sprmc
if not beam=="None":  catecut["beam"]  = beam
if not jbrsc=="None": catecut["jbrsc"] = jbrsc

# adjust
if EXTTRIG and DUMPROOT: ROOTEXT=True
if GRPTRIG and DUMPROOT: ROOTGRP=True

print ""
print "--- [OPTIONS] ----------------------"
print "    %s      "%(RUNTAG)
print "  (standard)"
print "    FORCE      = ", FORCE
print "    SUMMARY    = ", SUMMARY
print "    CALIB      = ", CALIB
print "    EXTTRIG    = ", EXTTRIG
print "    GRPTRIG    = ", GRPTRIG
print "    DELETE     = ", DELETE
print "    catecut    = ", catecut
print "    cut_pre    = ", cut_pre
print "    cut_post   = ", cut_post
print ""
print "  (ROOT)         "
print "    DUMPROOT   = ", DUMPROOT
print "    ROOTEXT    = ", ROOTEXT
print "    ROOTGRP    = ", ROOTGRP
print "------------------------------------"

npar = len(args)
seps=['_',',','-']
run_p=None
if (npar>=1):
    for sep in seps:
        if sep in str(args[0]): run_p = str(args[0]).split(sep)
    if run_p is None:
        run_p = str(args[0])
    print "analyzing runs: ",run_p
else:
    print "Error: specify run number of MuTES ", args
    sys.exit(0)

df = pd.read_csv(RUNINFO)
run_list = df.iloc[:,0].tolist()
noise_list = df.iloc[:,1].tolist()
ana_list = df.iloc[:,2].tolist()
exttrig_list = df.iloc[:,3].tolist()
grptrig_list = df.iloc[:,4].tolist()
cal_list = df.iloc[:,9].tolist()
if isinstance(run_p,list)==False:
    run_p = (run_p,)
    run_p = tuple(run_p)
run_p = [int(r) for r in run_p]
ind = run_list.index(run_p[0])
irn = int(noise_list[ind])
run_n="%04d"%irn
ana_target = ana_list[ind]
exttrig = exttrig_list[ind]
grptrig = grptrig_list[ind]

if ana_target=="":
    print "Error: ana is empty"
    sys.exit(0)
elif ana_list[ind]=="noise":
    print "Error: this is noise run, please select pulse run"
    sys.exit(0)

cal_run = None if str(cal_list[ind]) == "None" else int(cal_list[ind])
cal_noise_run = None
if cal_run is not None:
    calind = run_list.index(cal_run)
    cal_noise_run = int(noise_list[calind])
analist = [(run_p, int(run_n), ana_target, exttrig, grptrig, cal_run, cal_noise_run,BADCHS)]

for pulse_runnums, noise_runnum, target, extflag, grpflag, cal_runnum, cal_noisenum, badchan in analist:
    print "..... run, noise, target, exttrig, grptrig, cal_run, cal_noise = ", pulse_runnums, noise_runnum, target, extflag, grpflag, cal_runnum, cal_noisenum
    print "BADCHAN = ", badchan
    # ---------------------------------------------------------
    # when external trigger or group trigger is off, those flags are forced to be false. 
    if extflag == "off": 
    	orgEXTTRIG = EXTTRIG
        EXTTRIG = False
        ROOTEXT = False
        catecut["beam"] = "None"
    if grpflag == "off":
        orgGRPTRIG = GRPTRIG
        GRPTRIG = False
        ROOTGRP = False
        catecut["sprmc"] = "None"
        catecut["prime"] = "None"
    if extflag == "off" and grpflag == "off":
        catecut["jbrsc"] = "None"
#    if grouptrigmax: 
#        GRTINFO = "./csv/grptrig_singlethread_all.txt"
    # ---------------------------------------------------------
    mutes = MUTES.MUTES(pulse_runnums, noise_runnum, maxchans, cal_runnum, cal_noisenum,
                        badchan, DATADIR, DELETE, GRTINFO, COLUMN_INFO,
                        catecut=catecut, target=target,cut_pre=cut_pre, cut_post=cut_post)
    mutes.ana(forceNew=FORCE,summaryNew=SUMMARY,calibNew=CALIB,exttrigNew=EXTTRIG,grptrigNew=GRPTRIG)
    # ---------------------------------------------------------
    if DUMPROOT:
        print "\n [dump ROOT except for pulses]"
        mutes.dump_ROOT(EXTTRIG=ROOTEXT, GRTRIG=ROOTGRP)
    # ---------------------------------------------------------
    # when external trigger or group trigger is off, those flags are back to the input values 
    if extflag == "off":
        EXTTRIG = orgEXTTRIG
        catecut["beam"] = beam
    if grpflag == "off":
        GRPTRIG = orgGRPTRIG
        catecut["sprmc"] = sprmc
        catecut["prime"] = prime
    if extflag == "off" and grpflag == "off":
        catecut["jbrsc"] = jbrsc
    # adjust
    if EXTTRIG and DUMPROOT: ROOTEXT=True
    if GRPTRIG and DUMPROOT: ROOTGRP=True
    # ---------------------------------------------------------

