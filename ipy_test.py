#!/usr/bin/env python

""" run_mutes_single.py is a basic analysis tools for MUTES project. 

History: 
2018-07-04 ; ver 1.0; made from ana_galen_single.py, S.Yamada
2018-07-05 ; ver 1.1; update to dump ROOT and pulses 
2018-07-06 ; ver 1.2; a minor update
2019-01-31 ; ver 1.3; H.Tatsuno modified drastically, removed all plots functions and simplified
2019-04-26 ; ver 1.4; H.Tatsuno added multiple runs analysis
2019-05-09 ; ver 1.5; H.Tatsuno minor change

"""

__author__ =  'a bad boy HT'
__version__ = '1.5'

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

import mutes_ana as MUTES
MUTES= reload(MUTES) # for ipython to reload when it changed

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
#maxchans = 50
#print "maxchans is not full, this is a test mode, please do not use grp trigger"

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
if not jbrsc=="None": catecut["jbrsc"] = jbrsc
if not beam=="None":  catecut["beam"]  = beam
if not sprmc=="None": catecut["sprmc"] = sprmc
if beam=="None" and sprmc=="None": catecut=None

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

if isinstance(run_p,list)==False:
    run_p = (run_p,)
    run_p = tuple(run_p)

ind = run_list.index(run_p[0])
run_p = [int(r) for r in run_p]
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
    print "..... run, noise, target, exttrig, grptrig, cal_run, cal_noise = ", pulse_runnums, noise_runnum, target, extflag, grpflag, calibration_runnum, calibration_noisenum
    print "BADCHAN = ", badchan
    # when external trigger or group trigger is off, those flags are forced to be false. 
    if extflag == "off": 
    	orgEXTTRIG = EXTTRIG
        EXTTRIG = False
    if grpflag == "off":
        orgGRPTRIG = GRPTRIG
        GRPTRIG = False
#    if grouptrigmax: 
#        GRTINFO = "./csv/grptrig_singlethread_all.txt"
    # ---------------------------------------------------------
    mutes = MUTES.MUTES(pulse_runnums, noise_runnum, maxchans, calibration_runnum, calibration_noisenum,
                        badchan, DATADIR, DELETE, GRTINFO, COLUMN_INFO,
                        hdf5optname=hdf5optname, catecut=catecut, target=target,
                        cut_pre=cut_pre, cut_post=cut_post)
    mutes.ana(forceNew=FORCE,summaryNew=SUMMARY,calibNew=CALIB,exttrigNew=EXTTRIG,grptrigNew=GRPTRIG)
    # ---------------------------------------------------------
    #if DUMPROOT:
    #    print "\n [dump ROOT except for pulses]"
    #    mutes.dump_ROOT_2019(EXTTRIG=ROOTEXT, GRTRIG=ROOTGRP)
    # ---------------------------------------------------------
    # when external trigger or group trigger is off, those flags are back to the input values 
    if extflag == "off":    	EXTTRIG = orgEXTTRIG
    if grpflag == "off":    	GRPTRIG = orgGRPTRIG
    # ---------------------------------------------------------


# all data
data = mutes.data
chans=data.channel.keys()

## data set of first good channel
#ds1 = data.first_good_dataset
## data set of channel11 (if exists)
#ds11 = data.channel[11]
#
#ds = data.channel[99]
## good event array
#g=ds.good()
#
## filtered values
#filt_value = ds.p_filt_value
## with drift correction
#filt_value_dc = ds.p_filt_value_dc
## with drift correction + phase correction
#filt_value_phc = ds.p_filt_value_phc
## if calibration is succeeded
#energy        = ds.p_energy
## others
#pretrig_mean  = ds.p_pretrig_mean
#pretrig_rms   = ds.p_pretrig_rms
#pulse_average = ds.p_pulse_average
#peak_time     = ds.p_peak_time
#peak_value    = ds.p_peak_value
#rise_time     = ds.p_rise_time
#rowcount      = ds.p_rowcount
#timestamp     = ds.p_timestamp
#
## number of pulses
#npulse      = ds.nPulses
## number of samples (1024)
#nsamples    = ds.nSamples
## number of pre-samples (256)
#npresamples = ds.nPresamples
#
## pulse data (1024)
#pulse0 = ds.read_trace(0)
#pulse1 = ds.read_trace(1)
##.....
#pulse_last = ds.read_trace(npulse-1)
#
#
## how to select good event
## just extract indices with [g] g=ds.good()
#print "filt_value only good events", filt_value[g]
#print "energy only good events", energy[g]
#
## how to plot pulse
#plt.close('all')
#plt.ion()
#
#x=np.arange(nsamples)
#x=x-npresamples
#good_pulse_idx = np.where(g)[0]
#bad_pulse_idx = np.where(~g)[0]
#plt.figure()
#for i in range(5):
#    plt.plot(x,ds.read_trace(good_pulse_idx[i]))
#
#plt.figure()
#for i in range(5):
#    plt.plot(x,ds.read_trace(bad_pulse_idx[i]))
#
#plt.figure()
#plt.plot(timestamp[g],pretrig_mean[g])
#


# ------------------------------------------------------------

def localfit(ds,linename):
    '''
     NOTE parameters of MultiLorentzianComplexFitter
        param_meaning = {
            "resolution": 0,
            "peak_ph": 1,
            "dP_dE": 2,
            "amplitude": 3,
            "background": 4,
            "bg_slope": 5,
            "tail_frac": 6,
            "tail_length": 7
        }
    '''
    elo,ehi = mass.STANDARD_FEATURES[linename]-50,mass.STANDARD_FEATURES[linename]+50
    edges = np.arange(elo,ehi,1)
    counts, _ = np.histogram(ds.p_energy[ds.good()],edges)
    fitter = mass.getfitter(linename)
    fitter.fit(counts,edges,plot=False)
    params = fitter.last_fit_params[:]
    params[fitter.param_meaning["tail_frac"]]=0.25
    params[fitter.param_meaning["dP_dE"]]=1
    params[fitter.param_meaning["resolution"]]=6
    fitter.fit(counts,edges,params=params,hold=[2],vary_tail=True,plot=False)
    # hold[2]: fixing dP_dE
    #print "---------------------------"
    #print fitter.result_string()
    #print "---------------------------"
    return fitter


def decay_time_rough(ds,idx):
    y = ds.read_trace(idx)-ds.p_pretrig_mean[idx]
    peak_value = ds.p_peak_value[idx]
    peak_index = ds.p_peak_index[idx]
    dp = peak_value/np.exp(1)
    y2 = y[peak_index:]
    dis = np.where(np.abs(y2-dp)<50)[0]
    if len(dis==0):
        dis = np.where(np.abs(y2-dp)<100)[0]
    dt = np.mean(dis)*ds.timebase*1e6# micro sec
    return dt
    

# ------------------------------------------------------------

# how to fit
if ana_target=="Mn":
    linename="MnKAlpha"
elif ana_target=="Fe" or ana_target=="Co57":
    linename="FeKAlpha"

from matplotlib.backends.backend_pdf import PdfPages
import csv
import math

ncols = 8
nrows = 30
divx1, divy1 = 6, 5
divx2, divy2 = 2, 2

outdir="./output/"
fresolname=outdir+'run%d_resol.csv'%(run_p[0])

# energy resolution
#f = open(fresolname, 'w')
#writer = csv.writer(f, lineterminator='\n')
#res_list=[]
#pdfname=outdir+"run%d_fit_%s.pdf"%(run_p[0],linename)
#with PdfPages(pdfname) as pdf:
#    print "printing...", pdfname
#    for icol in xrange(1,ncols+1,1):
#        fig = plt.figure(figsize=(20,15))
#        for ich in xrange(1,nrows+1,1):
#            ch = (icol-1)*nrows*2+ich*2-1
#            if ch in chans:
#                ds = data.channel[ch]
#                fitter = localfit(ds,linename)
#                res = fitter.last_fit_params_dict["resolution"][0]
#                res_err = fitter.last_fit_params_dict["resolution"][1]
#                ch = ds.channum
#                print "%03d,%.2f,%.2f"%(ch,res,res_err)
#                if math.isnan(res):
#                    res=1e3
#                    res_err=1e3
#                res_list.append(["%03d"%ch,"%.2f"%res,"%.2f"%res_err])
#                ax = plt.subplot(divx1,divy1,ich)
#                fitter.plot(axis=ax,ph_units='eV')
#                ax.set_title("chan %d"%ch)
#        fig.tight_layout()
#        pdf.savefig()
#        plt.close()
#writer.writerows(res_list)
#f.close()
#print "%s is created."%pdfname
#print "%s is created."%fresolname

chs=[]
resols=[]
if os.path.isfile(fresolname):
    df = pd.read_csv(fresolname,header=None)
    chs    = df.iloc[:,0].tolist()
    resols = df.iloc[:,1].tolist()
    

pdfname=outdir+"run%d_pretrig_mean.pdf"%(run_p[0])
with PdfPages(pdfname) as pdf:
    print "printing...", pdfname
    for icol in xrange(1,ncols+1,1):
        fig = plt.figure(figsize=(20,15))
        for ich in xrange(1,nrows+1,1):
            ch = (icol-1)*nrows*2+ich*2-1
            if ch in chans:
                ds = data.channel[ch]
                ax = plt.subplot(divx1,divy1,ich)
                g=ds.good()
                ax.scatter(ds.p_timestamp[g],ds.p_pretrig_mean[g])
                ax.set_xlabel("timestamp")
                ax.set_ylabel("pretrig mean")
                ch = ds.channum
                print "%03d"%(ch)
                ax.set_title("chan %d"%ch)
                if ch in chs: ax.set_title("chan %d fwhm %.2f eV"%(ch,resols[chs.index(ch)]))
        fig.tight_layout()
        pdf.savefig()
        plt.close()



edges = np.arange(5000,15010,10)
pdfname=outdir+"run%d_filt_value.pdf"%(run_p[0])
with PdfPages(pdfname) as pdf:
    print "printing...", pdfname
    for icol in xrange(1,ncols+1,1):
        fig = plt.figure(figsize=(20,15))
        for ich in xrange(1,nrows+1,1):
            ch = (icol-1)*nrows*2+ich*2-1
            if ch in chans:
                ds = data.channel[ch]
                ax = plt.subplot(divx1,divy1,ich)
                counts, _ = np.histogram(ds.p_filt_value[ds.good()],edges)
                ax.step(edges[:-1],counts)
                ch = ds.channum
                print "%03d"%(ch)
                ax.set_title("chan %d filt_value"%ch)
                if ch in chs: ax.set_title("chan %d fwhm %.2f eV"%(ch,resols[chs.index(ch)]))
        fig.tight_layout()
        pdf.savefig()
        plt.close()


attr_phc='p_filt_value_phc'
attr_dc='p_filt_value_dc'
pdfname=outdir+"run%d_calibration_curve.pdf"%(run_p[0])
with PdfPages(pdfname) as pdf:
    print "printing...", pdfname
    for icol in xrange(1,ncols+1,1):
        fig = plt.figure(figsize=(20,15))
        for ich in xrange(1,nrows+1,1):
            ch = (icol-1)*nrows*2+ich*2-1
            if ch in chans:
                ds = data.channel[ch]
                ax = plt.subplot(divx1,divy1,ich)
                if ds.calibration.has_key(attr_phc): cal = ds.calibration[attr_phc]
                elif ds.calibration.has_key(attr_dc):cal = ds.calibration[attr_dc]
                else: continue
                cal.plot(axis=ax)
                ch = ds.channum
                print "%03d"%(ch)
                ax.set_title("chan %d energy calibration curve"%ch)
                if ch in chs: ax.set_title("chan %d fwhm %.2f eV"%(ch,resols[chs.index(ch)]))
        fig.tight_layout()
        pdf.savefig()
        plt.close()


x=np.arange(1024)
x=x-256
pdfname=outdir+"run%d_good_pulses.pdf"%(run_p[0])
with PdfPages(pdfname) as pdf:
    print "printing...", pdfname
    for icol in xrange(1,ncols+1,1):
        fig = plt.figure(figsize=(20,15))
        for ich in xrange(1,nrows+1,1):
            ch = (icol-1)*nrows*2+ich*2-1
            if ch in chans:
                ds = data.channel[ch]
                ax = plt.subplot(divx1,divy1,ich)
                good_pulse_idx = np.where(ds.good())[0]

                for i in range(2):
                    y = ds.read_trace(good_pulse_idx[i])-ds.p_pretrig_mean[good_pulse_idx[i]]
                    dt = decay_time_rough(ds,good_pulse_idx[i])
                    ax.plot(x*ds.timebase*1e3,y,label="decay time = %.0f us"%dt)
                    ax.set_xlabel("time (ms)")
                    ax.set_ylabel("fead back value (arb)")
                
                ch = ds.channum
                print "%03d"%(ch)
                ax.set_title("chan %d"%ch)
                if ch in chs: ax.set_title("chan %d fwhm %.2f eV"%(ch,resols[chs.index(ch)]))
                ax.legend()
        fig.tight_layout()
        pdf.savefig()
        plt.close()

