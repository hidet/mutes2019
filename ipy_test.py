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


# all data
data = mutes.data
chans=data.channel.keys()

# ------------------------------------------
## data set of first good channel
#ds1 = data.first_good_dataset
## data set of channel11 (if exists)
#ds11 = data.channel[11]
#
#ds = data.channel[23]
# good event array
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
#x=np.arange(1024)
#x=x-256
#good_pulse_idx = np.where(g)[0]
#bad_pulse_idx = np.where(~g)[0]
#plt.figure()
#for i in xrange(10):
#    y = ds.read_trace(good_pulse_idx[i])-ds.p_pretrig_mean[good_pulse_idx[i]]
#    plt.plot(x*ds.timebase*1e3,y)
#plt.title("chan %d"%ds.channum)
#plt.xlabel("Time (ms)")
#plt.ylabel("feed back (arb)")
#
#plt.figure()
#for i in xrange(30,40,1):
#    plt.plot(x*ds.timebase*1e3,ds.read_trace(good_pulse_idx[i]))
#plt.title("chan %d"%ds.channum)
#plt.xlabel("Time (ms)")
#plt.ylabel("feed back (arb)")
#
#plt.figure()
#for i in range(10):
#    plt.plot(x,ds.read_trace(bad_pulse_idx[i]))
#
#plt.figure()
#plt.plot(timestamp[g],pretrig_mean[g])
#
# ------------------------------------------

#sys.exit()


# ------------------------------------------------------------
# functions to run script as standalone
def localfit(ds,linename,category=None):
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
    if category is not None:
        good = ds.cuts.good(**category)
    else:
        good = ds.good()
    elo,ehi = mass.STANDARD_FEATURES[linename]-50,mass.STANDARD_FEATURES[linename]+50
    edges = np.arange(elo,ehi,1)
    counts, _ = np.histogram(ds.p_energy[good],edges)
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


def localfit_goodchs(data,linename,goodchs,category=None):
    elo,ehi = mass.STANDARD_FEATURES[linename]-50,mass.STANDARD_FEATURES[linename]+50
    edges = np.arange(elo,ehi,1)
    p_energy = []
    for ds in data:
        if not ds.channum in goodchs: continue
        if category is not None:
            good = ds.cuts.good(**category)
        else:
            good = ds.good()
        p_energy.extend(ds.p_energy[good])
    p_energy = np.asarray(p_energy, dtype=np.float32)
    counts, _ = np.histogram(p_energy,edges)
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



# ------------------------------------------------------------
from matplotlib.backends.backend_pdf import PdfPages
import csv
import math

# fit linename
if ana_target=="Mn":
    linename="MnKAlpha"
elif ana_target=="Fe" or ana_target=="Co57":
    linename="FeKAlpha"

# catecut adjust
if exttrig=="off" and grptrig=="off": catecut={}

# 8 columns, 30 rows, 240 chs
ncols = 8
nrows = 30
divx1, divy1 = 6, 5

outdir="./output/%s/"%(RUNTAG)
if os.path.isdir(outdir)==False: os.makedirs(outdir)
# ------------------------------------------------------------

# energy resolution
fresolname=outdir+'run%d_resol.csv'%(run_p[0])
if not os.path.isfile(fresolname):
    f = open(fresolname, 'w')
    writer = csv.writer(f, lineterminator='\n')
    res_list=[]
    pdfname=outdir+"run%d_fit_%s.pdf"%(run_p[0],linename)
    with PdfPages(pdfname) as pdf:
        print "printing...", pdfname
        for icol in xrange(1,ncols+1,1):
            fig = plt.figure(figsize=(20,15))
            for ich in xrange(1,nrows+1,1):
                ch = (icol-1)*nrows*2+ich*2-1
                if not ch in chans: continue
                ds = data.channel[ch]
                fitter = localfit(ds,linename,category=catecut)
                res = fitter.last_fit_params_dict["resolution"][0]
                res_err = fitter.last_fit_params_dict["resolution"][1]
                print "%03d,%.2f,%.2f"%(ch,res,res_err)
                if math.isnan(res) or res==0.:
                    res=1e3
                    res_err=1e3
                res_list.append(["%03d"%ch,"%.2f"%res,"%.2f"%res_err])
                ax = plt.subplot(divx1,divy1,ich)
                fitter.plot(axis=ax,ph_units='eV')
                ax.set_title("chan %d"%ch)
            fig.tight_layout()
            pdf.savefig()
            plt.close()
    writer.writerows(res_list)
    f.close()
    print "%s is created."%pdfname
    print "%s is created."%fresolname


chs=[]
resmax=10.# FWHM < 10 eV
goodchs=[]
resols=[]
if os.path.isfile(fresolname):
    df = pd.read_csv(fresolname,header=None)
    df.columns = ['ch','res','res_err']
    chs    = df.ch.tolist()
    resols = df.res.tolist()
    goodchs = df[df.res<resmax].ch.tolist()
    pdfname=outdir+"run%d_fit_%s_resolmap.pdf"%(run_p[0],linename)
    with PdfPages(pdfname) as pdf:
        ax=tesmap.plotDensityPlaneFromCSV(fresolname)
        pdf.savefig()
        plt.close()
    print "%s is created."%pdfname
    

pdfname=outdir+"run%d_fit_%s_%dpixels.pdf"%(run_p[0],linename,len(goodchs))
with PdfPages(pdfname) as pdf:
    fig = plt.figure(figsize=(8, 8))
    fitter = localfit_goodchs(data,linename,goodchs=goodchs,category=catecut)
    ref_mean = mass.STANDARD_FEATURES[linename]
    res = fitter.last_fit_params_dict["resolution"][0]
    res_err = fitter.last_fit_params_dict["resolution"][1]
    mean = fitter.last_fit_params_dict["peak_ph"][0]
    mean_err = fitter.last_fit_params_dict["peak_ph"][1]
    tailf = fitter.last_fit_params_dict["tail_frac"][0]
    tailf_err = fitter.last_fit_params_dict["tail_frac"][1]
    tailt = fitter.last_fit_params_dict["tail_length"][0]
    tailt_err = fitter.last_fit_params_dict["tail_length"][1]
    print "---- %s %d pixels----"%(linename,len(goodchs))
    print "ref mean:  %.4f         eV"%(ref_mean)
    print "fit mean:  %.4f +- %.4f eV"%(mean,mean_err)
    print "diff mean: %.4f +- %.4f eV"%(mean-ref_mean,mean_err)
    print "resol:     %.4f +- %.4f eV"%(res,res_err)
    print "fit tailf: %.4f +- %.4f   "%(tailf,tailf_err)
    print "fit tailt: %.4f +- %.4f eV"%(tailt,tailt_err)
    ax = plt.subplot(1,1,1)
    fitter.plot(axis=ax,ph_units='eV',label='full')
    ax.set_title("fit %s %d pixels (FWHM < %.1f eV)"%(linename,len(goodchs),resmax))
    fig.tight_layout()
    pdf.savefig()
    plt.close()
print "%s is created."%pdfname


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
                print "%d,"%(ch),
                sys.stdout.flush()
                ax.set_title("chan %d energy calibration curve"%ch)
                if ch in chs: ax.set_title("chan %d fwhm %.2f eV"%(ch,resols[chs.index(ch)]))
        fig.tight_layout()
        pdf.savefig()
        plt.close()
print "%s is created."%pdfname


        
pdfname=outdir+"run%d_pretrig_mean.pdf"%(run_p[0])
with PdfPages(pdfname) as pdf:
    print "printing...", pdfname
    for icol in xrange(1,ncols+1,1):
        fig = plt.figure(figsize=(20,15))
        for ich in xrange(1,nrows+1,1):
            ch = (icol-1)*nrows*2+ich*2-1
            if ch in chans:
                ds = data.channel[ch]
                g = ds.cuts.good(**catecut)# good event cuts
                ax = plt.subplot(divx1,divy1,ich)
                ax.scatter(ds.p_timestamp[g],ds.p_pretrig_mean[g])
                ax.set_xlabel("timestamp")
                ax.set_ylabel("pretrig mean")
                print "%d,"%(ch),
                sys.stdout.flush()
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
                g = ds.cuts.good(**catecut)# good event cuts
                ax = plt.subplot(divx1,divy1,ich)
                counts, _ = np.histogram(ds.p_filt_value[g],edges)
                ax.step(edges[:-1],counts)
                print "%d,"%(ch),
                sys.stdout.flush()
                ax.set_title("chan %d filt_value"%ch)
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
                g = ds.cuts.good(**catecut)# good event cuts
                ax = plt.subplot(divx1,divy1,ich)
                good_pulse_idx = np.where(g)[0]
                for i in range(2):
                    y = ds.read_trace(good_pulse_idx[i])-ds.p_pretrig_mean[good_pulse_idx[i]]
                    dt = decay_time_rough(ds,good_pulse_idx[i])
                    ax.plot(x*ds.timebase*1e3,y,label="decay time = %.0f us"%dt)
                    ax.set_xlabel("time (ms)")
                    ax.set_ylabel("fead back value (arb)")
                print "%d,"%(ch),
                sys.stdout.flush()
                ax.set_title("chan %d"%ch)
                if ch in chs: ax.set_title("chan %d fwhm %.2f eV"%(ch,resols[chs.index(ch)]))
                ax.legend()
        fig.tight_layout()
        pdf.savefig()
        plt.close()

