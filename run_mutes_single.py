#!/usr/bin/env python

""" run_mutes_single.py is a basic analysis tools for MUTES project. 

History: 
2018-07-04 ; ver 1.0; made from ana_galen_single.py, S.Yamada
2018-07-05 ; ver 1.1; update to dump ROOT and pulses 
2018-07-06 ; ver 1.2; a minor update
2019-01-31 ; ver 1.3; H.Tatsuno modified drastically, removed all plots functions and simplified
2019-04-26 ; ver 1.4; H.Tatsuno added multiple runs analysis

"""

__author__ =  'a bad boy HT'
__version__ = '1.4'

import matplotlib
matplotlib.use("Agg")
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

import mutes_ana as KHE
KHE= reload(KHE) # for ipython to reload when it changed

print "[START] " + __file__

ANADIR=os.environ.get("MUTESANADIR","")
DATADIR=os.environ.get("MUTESDATADIR","")
print "ANADIR  = ", ANADIR
print "DATADIR = ", DATADIR

if ANADIR == "" or DATADIR == "":
    print "[ERROR] Set MUTESANADIR and MUTESDATADIR"
    print "e.g., for bash users"
    sys.exit()
    
BADCHS = [3,9,39,77,83,85,111,337,367,375,423]# initially disconnected
BADCHS.extend([117,203,233])# bad channels
BADCHS.extend([5,177])# strange channels
BADCHS.extend([17])# ?? is this strange?
BADCHS.sort()

maxchans = 240
#maxchans = 10

if os.path.isdir(DATADIR)==False: 
    print "%s is missing"%DATADIR
    sys.exit()

RUNINFO="./csv/data_TMU_2019G.csv"
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


usage = u'(1) %prog 278 (basic analysis), (2) %prog 278'
version = __version__
parser = optparse.OptionParser(usage=usage, version=version)
# setting for options
# how to use
# -eg -REG : without forceNew, no summaryNew, no calibNew, with exttrig, grouptrig, and ROOT, ROOTEXT, ROOTGRP
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
parser.add_option('--hdf5optname',  dest='hdf5optname', action="store", type=str, help='add optional name for hdf5 (default=None)', default=None)

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
hdf5optname    = options.hdf5optname

catecut = {}
catecut["prime"] = "on"
if not beam=="None":
    catecut["beam"] = beam
if not sprmc=="None":
    catecut["sprmc"] = sprmc
if not jbrsc=="None":
    catecut["jbrsc"] = sprmc
if beam=="None" and sprmc=="None":
    if jbrsc=="None":
        catecut=None

print ""
print "--- [OPTIONS] ----------------------"
print "  (standard)"
print "    FORCE          = ", FORCE
print "    SUMMARY        = ", SUMMARY
print "    CALIB          = ", CALIB
print "    EXTTRIG        = ", EXTTRIG
print "    GRPTRIG        = ", GRPTRIG
print "    DELETE         = ", DELETE
print "    catecut        = ", catecut
print "    cut_pre        = ", cut_pre
print "    cut_post       = ", cut_post
print "    hdf5optname    = ", hdf5optname
print ""
print "  (ROOT)             "
print "    DUMPROOT       = ", DUMPROOT
print "    ROOTEXT        = ", ROOTEXT
print "    ROOTGRP        = ", ROOTGRP
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
if isinstance(run_p,list):
    ind = run_list.index(run_p[0])
    run_p = [int(r) for r in run_p]
else:
    ind = run_list.index(run_p)
    run_p = int(run_p)
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
    calind = run_list.index(str(cal_run))
    cal_noise_run = int(noise_list[calind])
analist = [(run_p, int(run_n), ana_target, exttrig, grptrig, cal_run, cal_noise_run,BADCHS)]

for pulse_runnums, noise_runnum, target, extflag, grpflag, calibration_runnum, calibration_noisenum, badchan in analist:
    print "..... run, noise, target, exttrig, grptrig, cal_run, = ", pulse_runnums, noise_runnum, target, extflag, grpflag, calibration_runnum, calibration_noisenum
    print "BADCHAN = ", badchan
    if extflag == "off": # when external trigger or group trigger is off, those flags are forced to be false. 
    	orgEXTTRIG = EXTTRIG
        EXTTRIG = False
    if grpflag == "off":
        orgGRPTRIG = GRPTRIG
        GRPTRIG = False
#    if grouptrigmax: 
#        GRTINFO = "./csv/grptrig_singlethread_all.txt"
    # ---------------------------------------------------------
    k = KHE.KHE(pulse_runnums, noise_runnum, maxchans, calibration_runnum, calibration_noisenum,
                badchan, DATADIR, DELETE, GRTINFO, COLUMN_INFO,
                hdf5optname=hdf5optname, catecut=catecut, target=target, cut_pre=cut_pre, cut_post=cut_post)
    k.anahide(forceNew=FORCE,summaryNew=SUMMARY,calibNew=CALIB,exttrigNew=EXTTRIG,grptrigNew=GRPTRIG)
    # ---------------------------------------------------------
    if DUMPROOT:
        print "\n [dump ROOT except for pulses]"
        k.dump_ROOT_2019(EXTTRIG=ROOTEXT, GRTRIG=ROOTGRP)
    # ---------------------------------------------------------
    if extflag == "off": # when external trigger or group trigger is off, those flags are back to the input values 
    	EXTTRIG = orgEXTTRIG
    if grpflag == "off":
    	GRPTRIG = orgGRPTRIG
    # ---------------------------------------------------------


# end
