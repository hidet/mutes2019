""" mutes_ext.py is a function package for MUTES E62 external trigger analysis 

History: 
2018-08-16 ; ver 1.0; branched from mutes_util.py, T.H.
2019-02-19 ; ver 1.1; modified by HT, especially the sign of adjust is opposite, removed round
2019-02-22 ; ver 1.2; modified by HT, previous modification was WRONG!! 
2019-05-09 ; ver 1.3; HT, spill was removfed for MUTES
2019-05-21 ; ver 1.4; HT, minor change

"""

__version__ = '1.4'

import mass
import numpy as np
import pylab as plt
import os
import h5py
import mutes_util as util
from collections import OrderedDict
import sys


def check_external_trigger_data(pulse_runnum,DATADIR):
    datadir = os.path.expanduser(DATADIR)
    dir_p = os.path.join(datadir,"run%04d"%pulse_runnum)
    return os.path.isfile("%s/run%04d_extern_trig.hdf5"%(dir_p,pulse_runnum))

def define_beam_timing(ds,external_trigger_rowcount,forceNew=False):
    if forceNew or not "beam" in ds.hdf5_group:
        p_row = np.asarray(ds.p_rowcount, dtype=np.int64) - util.GLOBAL_PT_OFFSET
        print "ch%d define beam timing ..."%(ds.channum)
        rows_after_last_external_trigger, rows_until_next_external_trigger\
            =mass.core.analysis_algorithms.nearest_arrivals(p_row[:], external_trigger_rowcount)
        rows_from_nearest_external_trigger = np.fmin(rows_after_last_external_trigger[:], rows_until_next_external_trigger[:])
        beamOn = rows_from_nearest_external_trigger*util.ROW_TIMEBASE<100e-06# within 100 us?
        ds.cuts.cut_categorical("beam", {"on": beamOn,"off": ~beamOn})
        h5_beam = ds.hdf5_group.require_dataset("beam",(ds.nPulses,),dtype=np.int32)
        h5_beam[:] = beamOn.astype(int)
    # becareful: the categorical cut means logical_and of good and beamOn
    try:
        print "loading beamflag hdf5 file ch%d..."%(ds.channum)
        setattr(ds,"p_beamflag",ds.hdf5_group["beam"][()])
    except:
        print "cannot load, no group of beam in hdf5 group"
    return

def calc_external_trigger_timing(ds,external_trigger_rowcount,forceNew=False):
    params=["rows_after_last_external_trigger_nrp",
            "rows_until_next_external_trigger_nrp",
            "rows_after_last_external_trigger_nrn",
            "rows_until_next_external_trigger_nrn"]
    if forceNew or not params[0] in ds.hdf5_group:
        print "ch%d external trigger analysis..."%(ds.channum)
        rows_after_last_external_trigger_nrp, rows_until_next_external_trigger_nrp\
            =mass.core.analysis_algorithms.nearest_arrivals(ds.p_rowp[:], external_trigger_rowcount)
        rows_after_last_external_trigger_nrn, rows_until_next_external_trigger_nrn\
            =mass.core.analysis_algorithms.nearest_arrivals(ds.p_rown[:], external_trigger_rowcount)
        g1nrp    = ds.hdf5_group.require_dataset(params[0],(ds.nPulses,), dtype=np.int64)
        g1nrp[:] = rows_after_last_external_trigger_nrp
        g2nrp    = ds.hdf5_group.require_dataset(params[1],(ds.nPulses,), dtype=np.int64)
        g2nrp[:] = rows_until_next_external_trigger_nrp
        g1nrn    = ds.hdf5_group.require_dataset(params[2],(ds.nPulses,), dtype=np.int64)
        g1nrn[:] = rows_after_last_external_trigger_nrn
        g2nrn    = ds.hdf5_group.require_dataset(params[3],(ds.nPulses,), dtype=np.int64)
        g2nrn[:] = rows_until_next_external_trigger_nrn
    try:
        print "loading ch%d hdf5 file..."%(ds.channum)
        for par in params:
            setattr(ds,par,ds.hdf5_group[par][()])
    except:
        print "cannot load, no groups in hdf5 group"
    return

def calc_spill_timing(ds,spill_start_rowcount,forceNew=False):
    params=["rows_after_last_spill_start",
            "rows_until_next_spill_start"]
    # for new_filter use ds.p_rowp, for old_filter use ds.p_rown
    if forceNew or not params[0] in ds.hdf5_group:
        print "ch%d spill timing analysis..."%(ds.channum)
        rows_after_last_spill_start, rows_until_next_spill_start\
            =mass.core.analysis_algorithms.nearest_arrivals(ds.p_rowp[:], spill_start_rowcount)
        g3       = ds.hdf5_group.require_dataset(params[0],(ds.nPulses,), dtype=np.int64)
        g3[:]    = rows_after_last_spill_start
        g4       = ds.hdf5_group.require_dataset(params[1],(ds.nPulses,), dtype=np.int64)
        g4[:]    = rows_until_next_spill_start
    try:
        print "loading ch%d hdf5 file..."%(ds.channum)
        for par in params:
            setattr(ds,par,ds.hdf5_group[par][()])
    except:
        print "cannot load, no groups in hdf5 group"
    return

#def calc_external_trigger_timing_from_dat(ds,external_trigger_rowcount,name='dat',forceNew=False):
#    if forceNew or not "rows_until_next_external_trigger_"+str(name) in ds.hdf5_group:
#        # ---------------------------------
#        print "ch%d external trigger analysis... %s"%(ds.channum,name)
#        rows_after_last_external_trigger, rows_until_next_external_trigger\
#            =mass.core.analysis_algorithms.nearest_arrivals(ds.p_nrp[:], external_trigger_rowcount)
#        g1 = ds.hdf5_group.require_dataset("rows_after_last_external_trigger_"+str(name),
#                                           (ds.nPulses,), dtype=np.int64)
#        g1[:] = rows_after_last_external_trigger
#        g2 = ds.hdf5_group.require_dataset("rows_until_next_external_trigger_"+str(name),
#                                           (ds.nPulses,), dtype=np.int64)
#        g2[:] = rows_until_next_external_trigger
#    else:
#        print "ch%d external trigger analysis skipped."%(ds.channum)
#    ds.rows_after_last_external_trigger_dat = ds.hdf5_group["rows_after_last_external_trigger_"+str(name)]
#    ds.rows_until_next_external_trigger_dat = ds.hdf5_group["rows_until_next_external_trigger_"+str(name)]
#    return
#
#def save_across_runs(runs,closeFigs=True,row=False,until=False):
#    dirname = "./output/across/{}".format(runs)
#    if not os.path.isdir(dirname):
#        os.makedirs(dirname)
#    for i in plt.get_fignums():
#        plt.figure(i)
#        if row:
#            if until: plt.savefig(os.path.join(dirname,"{}_row_until.png".format(i)))
#            else:     plt.savefig(os.path.join(dirname,"{}_row.png".format(i)))
#        else:
#            if until: plt.savefig(os.path.join(dirname,"{}_until.png".format(i)))
#            else:     plt.savefig(os.path.join(dirname,"{}.png".format(i)))
#        if closeFigs: plt.close()
#
#def save_across_runs_merged(runs,cut=False,closeFigs=True):
#    dirname = "./output/across/{}_{}".format(runs[0],runs[-1])
#    if not os.path.isdir(dirname):
#        os.makedirs(dirname)
#    for i in plt.get_fignums():
#        plt.figure(i)
#        if not cut:
#            plt.savefig(os.path.join(dirname,"{}.png".format(i)))
#        else:
#            plt.savefig(os.path.join(dirname,"{}_with_cut.png".format(i)))
#        if closeFigs: plt.close()        
#
#
#def hist2d_plot(X,Y,counts,shift,row,until,titlestr):
#    plt.figure()
#    if shift:
#        plt.pcolormesh(X-util.SHIFT_ROWS*util.ROW_TIMEBASE*1e6,Y,counts)
#    else:
#        plt.pcolormesh(X,Y,counts)
#    xl="hogehoge"
#    if row:
#        if not until: xl = "row after last external_trigger (row)"
#        else:         xl = "row until next external_trigger (row)"
#    else:
#        if not until: xl = "time after last external_trigger (us)"
#        else:         xl = "time until next external_trigger (us)"
#    plt.xlabel(xl)
#    plt.ylabel("energy (eV)")
#    plt.colorbar()
#    if shift:
#        plt.title("{}".format(titlestr))
#    else:
#        plt.title("{} SHIFTED".format(titlestr))
#
#
#def time_spectrum_plot(X,Y,counts,shift,row,until,titlestr):
#    plt.figure()
#    timecounts = counts.sum(axis=0)
#    timex = 0.5*(X[0,1:]+X[0,:-1])
#    if shift:
#        plt.plot(timex-util.SHIFT_ROWS*util.ROW_TIMEBASE*1e6, timecounts)
#    else:
#        plt.plot(timex, timecounts)
#    xl="hogehoge"
#    if row:
#        if not until: xl = "row after last external_trigger (row)"
#        else:         xl = "row until next external_trigger (row)"
#    else:
#        if not until: xl = "time after last external_trigger (us)"
#        else:         xl = "time until next external_trigger (us)"
#    plt.xlabel(xl)
#    plt.ylabel("counts per {:.2f} us bin".format(timex[1]-timex[0]))
#    if shift:
#        plt.title("{}\ntime spectrum".format(titlestr))
#    else:
#        plt.title("{} SHIFTED\ntime spectrum".format(titlestr))
#
#        
#def shiftXforplot(X,shift):
#    if shift:
#        Y=X-util.SHIFT_ROWS*util.ROW_TIMEBASE*1e6
#    else:
#        Y=X
#    return Y
#
## befor and after time differences are merged
#def hist2d_merged_plot(X,Y,counts,shift,titlestr):
#    plt.figure()
#    if shift:
#        plt.pcolormesh(X-util.SHIFT_ROWS*util.ROW_TIMEBASE*1e6,Y,counts)
#    else:
#        plt.pcolormesh(X,Y,counts)
#    xl="row count difference (1row=240ns)"
#    plt.xlabel(xl)
#    plt.ylabel("energy (eV)")
#    plt.colorbar()
#    if not shift:
#        plt.title("{}".format(titlestr))
#    else:
#        plt.title("{} SHIFTED".format(titlestr))
#
#def hist2d_merged(ds,energy_edges,time_edges,shift,category={}):
#    if not shift:
#        x1=ds.rows_after_last_external_trigger[ds.good(**category)]
#        x2=-ds.rows_until_next_external_trigger[ds.good(**category)]
#    else:
#        x1=ds.rows_after_last_external_trigger[ds.good(**category)]
#        x2=-ds.rows_until_next_external_trigger[ds.good(**category)]
#    y=ds.p_energy[ds.good(**category)]
#    counts1, _, _ = np.histogram2d(x1,y,(time_edges,energy_edges))
#    counts2, _, _ = np.histogram2d(x2,y,(time_edges,energy_edges))
#    counts = counts1 + counts2
#    counts = counts.T
#    X, Y = np.meshgrid(time_edges, energy_edges)
#    return X,Y,counts
#
#def hist2d_data_merged(data,elo,ehi,ebin,tlo,thi,tbin,shift,plot,category={}):
#    time_edges = np.arange(tlo,thi,tbin)
#    energy_edges = np.arange(elo,ehi,ebin)
#    counts = np.zeros(shape=(len(energy_edges)-1,len(time_edges)-1))
#    for ds in data:
#        X,Y,c = hist2d_merged(ds,energy_edges,time_edges,shift,category)
#        counts+=c
#    if plot:
#        hist2d_merged_plot(X,Y,counts,shift,os.path.split(ds.filename)[-1][:7]+", {} good chans".format(data.num_good_channels))
#    return X,Y,counts
#
#def hist2d_merged_cut(ds,energy_edges,time_edges,shift,cutarray):
#    if not shift:
#        x1=ds.rows_after_last_external_trigger[cutarray]
#        x2=-ds.rows_until_next_external_trigger[cutarray]
#    else:
#        x1=ds.rows_after_last_external_trigger[cutarray]
#        x2=-ds.rows_until_next_external_trigger[cutarray]
#    y=ds.p_energy[cutarray]
#    counts1, _, _ = np.histogram2d(x1,y,(time_edges,energy_edges))
#    counts2, _, _ = np.histogram2d(x2,y,(time_edges,energy_edges))
#    counts = counts1 + counts2
#    counts = counts.T
#    X, Y = np.meshgrid(time_edges, energy_edges)
#    return X,Y,counts
#
#def hist2d_data_merged_cut(data,elo,ehi,ebin,tlo,thi,tbin,shift,plot,cutarray):
#    time_edges = np.arange(tlo,thi,tbin)
#    energy_edges = np.arange(elo,ehi,ebin)
#    counts = np.zeros(shape=(len(energy_edges)-1,len(time_edges)-1))
#    for onecut,ds in enumerate(data):
#        X,Y,c = hist2d_merged_cut(ds,energy_edges,time_edges,shift,cutarray[onecut])
#        counts+=c
#    if plot:
#        hist2d_merged_plot(X,Y,counts,shift,os.path.split(ds.filename)[-1][:7]+", {} good chans".format(data.num_good_channels))
#    return X,Y,counts
#
#def time_spectrum_merged_plot(X,Y,counts,shift,titlestr):
#    plt.figure()
#    timecounts = counts.sum(axis=0)
#    timex = 0.5*(X[0,1:]+X[0,:-1])
#    if shift:
#        plt.plot(timex-util.SHIFT_ROWS*util.ROW_TIMEBASE*1e6, timecounts)
#    else:
#        plt.plot(timex, timecounts)
#    xl="row count difference (1row=240ns)"
#    plt.xlabel(xl)
#    plt.ylabel("counts per {:.2f} row bin".format(timex[1]-timex[0]))
#    if shift:
#        plt.title("{}\ntime spectrum".format(titlestr))
#    else:
#        plt.title("{} SHIFTED\ntime spectrum".format(titlestr))
#
#def hist2d_across_runs(pruns,nruns,closeFigs,add):
#    d = OrderedDict()
#    for i,j in zip(pruns,nruns):
#        dirname = "./output/run{:04d}_n{:04d}".format(i,j)
#        print "reading %s"%dirname
#        h5name = os.path.join(dirname,"2dhist" + add + ".hdf5") 
#        if os.path.isfile(h5name):
#            d[i]=h5py.File(h5name,"r")
#        else:
#            raise Exception("{} does not exist".format(i))
#    suffixes = [u'_CrKAlpha_shift=False',
#     u'_CrKAlpha_shift=True',
#     u'_FeKAlpha_shift=False',
#     u'_FeKAlpha_shift=True',
#     u'_He3_shift=False',
#      u'_He3_shift=True',
#      u'_He4_shift=False',
#      u'_He4_shift=True']
#
#    outd = OrderedDict()
#    for suffix in suffixes:
#        h5 = d.values()[0]
#        Xall = h5["X{}".format(suffix)].value
#        Yall = h5["Y{}".format(suffix)].value
#        counts = h5["counts{}".format(suffix)].value
#        counts_sum = np.zeros_like(counts)
#
#        for (i,h5) in d.items():
#            X = h5["X{}".format(suffix)].value
#            Y = h5["Y{}".format(suffix)].value
#            counts = h5["counts{}".format(suffix)].value
#            if not np.array_equal(X,Xall):
#                raise Exception("X doesnt match for {}".format(suffix))
#            counts_sum+=counts
#        outd[suffix]=X,Y,counts_sum
#
#    for suffix,(X,Y,counts_sum) in outd.items():
#        shift = "True" in suffix
#        linename = suffix[1:8]
#        hist2d_plot(X,Y,counts_sum,shift,row,until,"runs {}\n{}".format(d.keys(),linename))
#        time_spectrum_plot(X,Y,counts_sum,shift,row,until,"runs {}\n{}".format(d.keys(),linename))
#    save_across_runs(pruns,closeFigs,row,until)
#    return outd
#
#def hist2d_across_runs_merged(pruns,nruns,cut,closeFigs):
#    d = OrderedDict()
#    for i,j in zip(pruns,nruns):
#        dirname = "./output/run{:04d}_n{:04d}".format(i,j)
#        print "reading %s"%dirname
#        if not cut:
#            h5name = os.path.join(dirname,"2dhist_row_merged.hdf5")
#        else:
#            h5name = os.path.join(dirname,"2dhist_row_merged_with_cut.hdf5")
#        if os.path.isfile(h5name):
#            d[i]=h5py.File(h5name,"r")
#        else:
#            raise Exception("{} does not exist".format(i))
#    suffixes = [u'_CrKAlpha_shift=False',
#     u'_CrKAlpha_shift=True',
#     u'_FeKAlpha_shift=False',
#     u'_FeKAlpha_shift=True',
#     u'_He3_shift=False',
#      u'_He3_shift=True',
#      u'_He4_shift=False',
#      u'_He4_shift=True']
#
#    outd = OrderedDict()
#    for suffix in suffixes:
#        h5 = d.values()[0]
#        Xall = h5["X{}".format(suffix)].value
#        Yall = h5["Y{}".format(suffix)].value
#        counts = h5["counts{}".format(suffix)].value
#        counts_sum = np.zeros_like(counts)
#
#        for (i,h5) in d.items():
#            X = h5["X{}".format(suffix)].value
#            Y = h5["Y{}".format(suffix)].value
#            counts = h5["counts{}".format(suffix)].value
#            if not np.array_equal(X,Xall):
#                raise Exception("X doesnt match for {}".format(suffix))
#            counts_sum+=counts
#        outd[suffix]=X,Y,counts_sum
#
#    for suffix,(X,Y,counts_sum) in outd.items():
#        shift = "True" in suffix
#        linename = suffix[1:8]
#        hist2d_merged_plot(X,Y,counts_sum,shift,"runs {}\n{}".format(d.keys(),linename))
#        time_spectrum_merged_plot(X,Y,counts_sum,shift,"runs {}\n{}".format(d.keys(),linename))
#    save_across_runs_merged(pruns,cut,closeFigs)
#    return outd
#
#def get_across_XYC(pruns,nruns,closeFigs,add):
#    d = OrderedDict()
#    for i,j in zip(pruns,nruns):
#        dirname = "./output/run{:04d}_n{:04d}".format(i,j)
#        print "reading %s"%dirname
#        h5name = os.path.join(dirname,"2dhist" + add + ".hdf5") 
#        if os.path.isfile(h5name):
#            d[i]=h5py.File(h5name,"r")
#        else:
#            raise Exception("{} does not exist".format(i))
#    suffixes =[
#        u'_FeKAlpha_shift=False',
#        u'_He3_shift=False',
#        u'_He4_shift=False'
#    ]
#    outd = OrderedDict()
#    for suffix in suffixes:
#        h5 = d.values()[0]
#        Xall = h5["X{}".format(suffix)].value
#        Yall = h5["Y{}".format(suffix)].value
#        counts = h5["counts{}".format(suffix)].value
#        counts_sum = np.zeros_like(counts)
#        for (i,h5) in d.items():
#            X = h5["X{}".format(suffix)].value
#            Y = h5["Y{}".format(suffix)].value
#            counts = h5["counts{}".format(suffix)].value
#            if not np.array_equal(X,Xall):
#                raise Exception("X doesnt match for {}".format(suffix))
#            counts_sum+=counts
#        outd[suffix]=X,Y,counts_sum
#    return outd
#
#def get_across_XYC_merged(pruns,nruns,closeFigs):
#    return get_across_XYC(pruns,nruns,closeFigs,"_row_merged")
#
#def get_across_XYC_merged_with_cut(pruns,nruns,closeFigs):
#    return get_across_XYC(pruns,nruns,closeFigs,"_row_merged_with_cut")
#
#def find_runs():
#    runs=[]
#    for i in range(10000):
#        dirname = "./output/run{:04d}".format(i)
#        h5name = os.path.join(dirname,"2dhist.hdf5")
#        if os.path.isfile(h5name):
#            runs.append(i)
#    return runs
#
#def plot_and_save_hist2d_merged_with_cut(data,fname,cutarray,ebin=1.0,tlo=0,thi=100,tbin=1):
#    print "%s is creating..."%fname
#    with h5py.File(fname,"w") as h5:
#        for shift in [True,False]:
#            for linename in ["FeKAlpha","CrKAlpha","He3","He4"]:
#                elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
#                X,Y,counts =hist2d_data_merged_cut(data,elo,ehi,ebin,tlo,thi,tbin,
#                                                       shift=shift,plot=False,cutarray=cutarray)
#                h5["X_{}_shift={}".format(linename,shift)]=X
#                h5["Y_{}_shift={}".format(linename,shift)]=Y
#                h5["counts_{}_shift={}".format(linename,shift)]=counts
#                
#        elo,ehi = np.array([2000,10000])
#        X,Y,counts =hist2d_data_merged_cut(data,elo,ehi,ebin,tlo,thi,tbin,
#                                               shift=shift,plot=False,cutarray=cutarray)
#        h5["X_all_shift={}".format(shift)]=X
#        h5["Y_all_shift={}".format(shift)]=Y
#        h5["counts_all_shift={}".format(shift)]=counts
#    print "%s is created."%fname
#
#def plot_and_save_hist2d_merged(data,fname,ebin=1.0,tlo=0,thi=100,tbin=1):
#    print "%s is creating..."%fname
#    with h5py.File(fname,"w") as h5:
#        for shift in [True,False]:
#            for linename in ["FeKAlpha","CrKAlpha","He3","He4"]:
#                elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
#                X,Y,counts =hist2d_data_merged(data,elo,ehi,ebin,tlo,thi,tbin,
#                                               shift=shift,plot=True)
#                h5["X_{}_shift={}".format(linename,shift)]=X
#                h5["Y_{}_shift={}".format(linename,shift)]=Y
#                h5["counts_{}_shift={}".format(linename,shift)]=counts
#                elo,ehi = np.array([0,15000])
#        X,Y,counts =hist2d_data_merged(data,elo,ehi,ebin,tlo,thi,tbin,
#                                       shift=shift,plot=True)
#        h5["X_all_shift={}".format(shift)]=X
#        h5["Y_all_shift={}".format(shift)]=Y
#        h5["counts_all_shift={}".format(shift)]=counts
#    print "%s is created."%fname


class Spill():
    MIN_SPILLS = 50
    #def __init__(self, ds, threshold=0.4):
    # 2018/06/20 Hide: changed to use self extrn values
    def __init__(self, ds, extrn_trig_rowc=None, extrn_trig_rowc_as_sec=None, threshold=0.4):
        self.ds = ds
        self.runstr = os.path.split(ds.filename)[-1][:7]
        if extrn_trig_rowc == None:
            try:
                self.external_trigger_rowcount = np.asarray(ds.external_trigger_rowcount[:], dtype=np.int64)
                self.external_trigger_rowcount_as_seconds = np.array(ds.external_trigger_rowcount_as_seconds[:])
            except:
                print "no spill, no external trigger data in ds"
                return None
        else:
            self.external_trigger_rowcount = extrn_trig_rowc
            self.external_trigger_rowcount_as_seconds = extrn_trig_rowc_as_sec
        self.spill_starts(threshold)
        if not self.have_enough_spills():
            print "number of spills is too small"
            return None
        
        self.spill_start_s = self.external_trigger_rowcount_as_seconds[np.array(self.spill_start_inds)]
        self.spill_start_rowcount = self.external_trigger_rowcount[:][np.array(self.spill_start_inds)]
        self.spill_end_s = self.external_trigger_rowcount_as_seconds[np.array(self.spill_end_inds)]
        self.spill_end_rowcount = self.external_trigger_rowcount[:][np.array(self.spill_end_inds)]      
        self.spill_period_s = np.median(np.diff(self.spill_start_s))
        self.typical_external_trigger_spacing = np.median(np.diff(self.external_trigger_rowcount_as_seconds))

    def have_enough_spills(self):
        return len(self.spill_start_inds)>self.MIN_SPILLS

    def spill_starts(self,threshold):
        # find spill structure t0
        # 1. high intensity beam -> lots of external triggers
        # 2. low intensity beam -> few external triggers
        state = 1
        spill_start_inds = []
        spill_end_inds = []
        for i,d in enumerate(np.diff(self.external_trigger_rowcount_as_seconds)):
            if state == 0:
                if d<threshold:
                    spill_start_inds.append(i)
                    state=1
            elif state == 1:
                if d>threshold:
                    spill_end_inds.append(i)
                    state = 0
        self.spill_start_inds = np.array(spill_start_inds[:len(spill_end_inds)-1],dtype=np.int64)
        self.spill_end_inds = np.array(spill_end_inds[1:],dtype=np.int64) # make sure there are the same number of starts and ends, and and is always after start

    def plot_spill_status(self):
        t = self.external_trigger_rowcount_as_seconds
        # lets plot the first 10 spills
        n_spill = 10
        n = self.spill_end_inds[n_spill]
        plt.figure()
        plt.plot(t[:n],t[:n]>0,".",label="external trigger")
        plt.plot([t[i] for i in self.spill_start_inds[:n_spill]], np.ones(n_spill),"r.",label="start")
        plt.plot([t[i] for i in self.spill_end_inds[:n_spill]], np.ones(n_spill),"m.",label="end")
        plt.xlabel("time (s)")
        plt.ylabel("spill status")
        plt.ylim(0,1.01)
        plt.legend(loc="best")
        plt.title("{}, chan {:d}".format(self.runstr,self.ds.channum))

    def plot_spill_start_differences(self):
        plt.figure()
        plt.plot(self.spill_start_s[1:],np.diff(self.spill_start_s),".")
        plt.xlabel("time (s)")
        plt.ylabel("duration between spill starts (s)")
        plt.title("{}, chan {:d}".format(self.runstr,self.ds.channum))

    def plot_spill_intensity(self):
        t = self.external_trigger_rowcount_as_seconds
        plt.figure()
        plt.plot(t[self.spill_start_inds[:-1]], np.diff(self.spill_start_inds),".")
        plt.xlabel("time (s)")
        plt.ylabel("external triggers during spill")
        plt.title("{}, chan {:d}".format(self.runstr,self.ds.channum))

    def ptmean_plot(self,spill_bin_s=0.01):
        spill_edges = np.arange(0,self.spill_period_s+spill_bin_s,spill_bin_s)
        spill_midpoints = 0.5*(spill_edges[1:]+spill_edges[:-1])
        self.rows_ext_after_last_spill_start, self.rows_ext_until_next_spill_start\
            =mass.core.analysis_algorithms.nearest_arrivals(self.external_trigger_rowcount,
                                                            np.array(self.spill_start_rowcount,dtype=np.int64))
        counts, _ = np.histogram(self.rows_ext_after_last_spill_start*util.ROW_TIMEBASE,spill_edges )
        plt.figure()
        inds = self.ds.rows_after_last_spill_start<np.iinfo(np.int64).max/10
        plt.plot(self.ds.rows_after_last_spill_start[inds]*util.ROW_TIMEBASE,self.ds.p_pretrig_mean[inds],".",label="pretrig mean")
        plt.xlim(0,5.2)
        a,b=np.percentile(self.ds.p_pretrig_mean,[10,90])
        plt.ylim(a-50,b+50)
        counts_plot = counts*1.0*(b-a)/np.amax(counts)+a
        plt.plot(spill_midpoints, counts_plot,label="beam current (arb)")
        plt.title("{}, chan {:d}".format(self.runstr,self.ds.channum))
        plt.legend(loc="best")
        plt.xlabel("time after spill start (s)")

    def plot_all(self,spill_bin_s=0.01):
        if self.have_enough_spills():
            self.plot_spill_start_differences()
            self.plot_spill_intensity()
            self.ptmean_plot(spill_bin_s=0.01)
            self.plot_spill_status()
        else:
            plt.figure()
            plt.title("only {} spills found, can't make spill plots".format(len(self.spill_start_inds)))

