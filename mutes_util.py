""" mutes_util.py is a function package for MUTES E62 analysis

History: 
2018-07-04 ; ver 1.0; created from khe.py, E62 run 
2018-07-05 ; ver 1.1; minor update by S.Y. 
2018-07-06 ; ver 1.2; minor update by S.Y. 

"""

__version__ = '1.2'

import mass
import numpy as np
import pandas as pd
import pylab as plt
params = {'xtick.labelsize': 10, # x ticks
          'ytick.labelsize': 10, # y ticks
          'legend.fontsize': 8
                    }
plt.rcParams['font.family'] = 'serif' # 
plt.rcParams.update(params)
import os
import h5py
from collections import OrderedDict
import sys

SHIFT_ROWS = 100
ROW_TIMEBASE = 7.18e-06/30.0
NUM_ROWS = 30

# global pulse timing offset
# (nSamples - npreSamples) * number of rows
GLOBAL_PT_OFFSET = (1024-256)*NUM_ROWS

hpht_phc="hpht_phc"# hname  

def get_usechans(pulse_runnum, noise_runnum, maxchans, badchans, DATADIR):
    datadir = os.path.expanduser(DATADIR)
    dir_p = os.path.join(datadir,"run%04d"%pulse_runnum)
    dir_n = os.path.join(datadir,"run%04d"%noise_runnum)
    avail_chans = mass.ljh_util.ljh_get_channels_both(dir_p,dir_n)
    chans = avail_chans[:min(len(avail_chans),maxchans)]
    for ch in badchans:
        if ch in chans: chans.remove(ch)
    return chans

def get_file_lists(pulse_runnum, noise_runnum, maxchans, badchans, DATADIR):
    datadir = os.path.expanduser(DATADIR)
    dir_p = os.path.join(datadir,"run%04d"%pulse_runnum)
    dir_n = os.path.join(datadir,"run%04d"%noise_runnum)
    avail_chans = mass.ljh_util.ljh_get_channels_both(dir_p,dir_n)
    chans = avail_chans[:min(len(avail_chans),maxchans)]
    for ch in badchans:
        if ch in chans: chans.remove(ch)
    pulse_files = [f+".ljh" for f in mass.ljh_util.ljh_chan_names(dir_p,chans)]
    noise_files = [f+".ljh" for f in mass.ljh_util.ljh_chan_names(dir_n,chans)]
    return pulse_files, noise_files

def ljh_fnames_ch(fnames, chan):
    basenames = [mass.ljh_util.ljh_basename_channum(fname)[0] for fname in fnames]
    return ["%s_chan%d"%(basename,chan) for basename in basenames]

def get_multiple_file_lists(pulse_runnums, noise_runnum, maxchans, badchans, DATADIR):
    datadir = os.path.expanduser(DATADIR)
    dir_ps = [os.path.join(datadir,"run%04d"%pulse_runnum) for pulse_runnum in pulse_runnums]
    dir_n = os.path.join(datadir,"run%04d"%noise_runnum)
    avail_chans = mass.ljh_util.ljh_get_channels_both(dir_ps[0],dir_n)
    chans = avail_chans[:min(len(avail_chans),maxchans)]
    for ch in badchans:
        if ch in chans: chans.remove(ch)
    pulse_files = [[f+".ljh" for f in ljh_fnames_ch(dir_ps,chan)] for chan in chans]
    noise_files = [f+".ljh" for f in mass.ljh_util.ljh_chan_names(dir_n,chans)]
    return pulse_files, noise_files

def generate_hdf5_filename(rawname, add="_mass"):
    """Generate the appropriate HDF5 filename based on a file's LJH name.
    Takes /path/to/data_chan33.ljh --> /path/to/data_mass.hdf5
    """
    import re
    fparts = re.split(r"_chan\d+", rawname)
    prefix_path = fparts[0]
    return prefix_path + add + ".hdf5"

def generate_root_filename(rawname, add="_mass"):
    """Generate the appropriate root filename based on a file's LJH name.
    Takes /path/to/data_chan33.ljh --> /path/to/data_mass.root
    """
    import re
    fparts = re.split(r"_chan\d+", rawname)
    prefix_path = fparts[0]
    return prefix_path + add + ".root"

def generate_user_root_filename(dname, add="_mass"):
    """generate root file name at the user specified directory
    """
    return dname + add + ".root"

def init_row_timebase(ds):
    global NUM_ROWS
    NUM_ROWS = ds.pulse_records.datafile.number_of_rows
    global ROW_TIMEBASE
    ROW_TIMEBASE = ds.timebase/NUM_ROWS
    global GLOBAL_PT_OFFSET
    GLOBAL_PT_OFFSET = (ds.nSamples - ds.nPresamples)*NUM_ROWS
    #print "GLOBAL PULSE TIME OFFSET = %d rows"%(GLOBAL_PT_OFFSET)

#def get_beam_runs(tesrun):
#    fname='csv/beam_tes_runs.csv'
#    if os.path.isfile(fname):
#        data = pd.read_table(fname,sep=' ',names=['tes','beam'])
#        beamruns=data.query('tes==%d'%tesrun)['beam'].values
#        return beamruns
#    else:
#        print "ERROR : couldn't open ", fname
#        sys.exit()
#
#def get_beam_clocks(prefix,tesrun):
#    beamruns=get_beam_runs(tesrun)
#    CLOCKDATADIR=os.environ.get("CLOCKDATADIR","")
#    clocks = np.asarray([],dtype=np.int64)
#    for run in beamruns:
#        filename="%s/%s_run%05d.dat" % (CLOCKDATADIR, prefix,run)
#        print 'beam clock file name: ' + filename
#        if os.path.isfile(filename):
#            tmp=np.asarray(np.loadtxt(filename),dtype=np.int64)
#            clocks=np.append(clocks,tmp)
#        else:
#            print "ERROR : couldn't open ", filename
##            sys.exit()
#    return clocks

# to get group triger map
def get_grtdict(fname):
    if os.path.isfile(fname):
        fin = open(fname,'r')
        "in readcsv, open ", fname
        gdict = {}
        for i, oneline in enumerate(fin):
            list = oneline.strip().split(',')
            list = [int(a) for a in list]
            gdict[list[0]] = list[1:]
        fin.close()
        return gdict
    else:
        print "ERROR : couldn't open ", fname
        sys.exit()


def get_column_info(fname):
    if os.path.isfile(fname):
        fin = open(fname,'r')
        "in readcsv, open ", fname
        col_dict = {}
        for i, oneline in enumerate(fin):
            if i == 0: continue
            list = oneline.strip().split(',')
            list = [a for a in list]
            col_dict[list[5]] = list[4]
        fin.close()
        return col_dict
    else:
        print "ERROR : couldn't open ", fname
        sys.exit()


def output_mux_neighbor_channels(colinfo):
    if os.path.isfile(colinfo) is False:
        print "ERROR : couldn't open ", colinfo
        sys.exit()
    import csv
    f = open('./csv/mux_neighbor.csv', 'w')
    writer = csv.writer(f, lineterminator='\n')
    
    df = pd.read_csv(colinfo)
    df = df.rename(columns={'Unnamed: 0': 'bay'})
    df = df.rename(columns={'Unnamed: 1': 'bay number'})
    baynames=['AX','AY','BX','BY','CX','CY','DX','DY']

    for bayname in baynames:
        dfbay = df[df['bay'].isin([bayname])]
        chs = dfbay['matter'].values
        rs = dfbay['RS'].values
        for i,ch in enumerate(chs):
            neichs = dfbay.query('RS==%d-1 or RS==%d+1'%(rs[i],rs[i]))['matter'].values
            outa=np.insert(neichs,0,ch)
            writer.writerow(outa.tolist())
    f.close()



def output_grpt_wo_mux_neighbor_channels(fgrp,fnei):
    if os.path.isfile(fgrp) is False or os.path.isfile(fnei) is False:
        print "ERROR : couldn't open ", fgrp, fnei
        sys.exit()
    import csv
    f = open('./csv/grptrig_wo_neighbor.csv', 'w')
    writer = csv.writer(f, lineterminator='\n')

    dfgrp = pd.read_csv(fgrp,header=None, index_col=0)
    dfnei = pd.read_csv(fnei,header=None, index_col=0)
    chs = dfgrp.index.values
    for i,ch in enumerate(chs):
        grplist=dfgrp.query('index==%d'%(ch)).values[0].tolist()
        neilist=dfnei.query('index==%d'%(ch)).values[0].tolist()
        grp_chs_new = [int(x) for x in grplist if not int(x) in neilist]
        grp_chs_new.insert(0,ch)
        writer.writerow(grp_chs_new)
    f.close()
    

    
def linefitds(ds, linename):
    elo,ehi = mass.STANDARD_FEATURES[linename]-50,mass.STANDARD_FEATURES[linename]+50
    edges = np.arange(elo,ehi,1)
    counts, _ = np.histogram(ds.p_energy[ds.good()],edges)
    fitter = getattr(mass,linename+"Fitter")()
    fitter.fit(counts, edges,plot=False,)
    return fitter

def linefit(counts, edges, linename):
    fitter = getattr(mass,linename+"Fitter")()
#    fitter = getattr(mass.fluorescence_lines,linename)()
    fitter.fit(counts, edges, plot=False)
    params = fitter.last_fit_params[:]
    params[fitter.param_meaning["tail_frac"]]=0.25
    params[fitter.param_meaning["dP_dE"]]=1
    params[fitter.param_meaning["resolution"]]=6
    # fit while holding dP_dE fixed
    fitter.fit(counts, edges,params,hold=[2],vary_tail=True, plot=False)
    return fitter

def subsample_arrival_plot(ds):
    # make a function to plot average pulse vs true_phase
    phase_edges = np.arange(-0.5,0.6,0.1)
    tavgs = OrderedDict()
    for i in range(len(phase_edges)-1):
        lo,hi = phase_edges[i],phase_edges[i+1]

        inds = np.nonzero(np.logical_and(ds.p_filt_phase[:]>lo,ds.p_filt_phase[:]<hi))[0]
        tavg = np.zeros(ds.nSamples-1)
        for i in inds:
            t=ds.read_trace(i)*1.0-ds.p_pretrig_mean[i]
            if ds.p_shift1[i]:
                tavg+=t[:-1]
            else:
                tavg+=t[1:]
        if len(inds)>5:
            tavg/=len(inds)
            tavgs[lo,hi]=tavg,len(inds)
    plt.figure()
    for (lo,hi),(tavg,npulses) in tavgs.items():
        print("{:.2f} {:.2f}".format(np.mean([lo,hi]),np.dot(tavg,ds.filter.filt_aterms[0])))

        plt.plot(tavg, label="{:.2f}, {:d}".format(np.mean([lo,hi]),npulses))
    plt.legend(title="p_filt_phase,npulses")
    plt.xlim(ds.nPresamples, ds.nPresamples+4)
    plt.ylim(-1,np.amax([tavg[ds.nPresamples+4] for tavg,npulses in tavgs.values()]))
    plt.grid(True)
    plt.xlabel("Sample Number")
    plt.ylabel("fb signal")
    runstr = os.path.split(ds.filename)[-1][:7]
    plt.title("{}, chan {:d}\navg of pulse[:-1] if shift1, pulse[1:] if not shift1".format(runstr,ds.channum))


    # make a function to plot average pulse vs true_phase
    phase_plus = -(ds.p_shift1[:]+ds.p_filt_phase[:])
    phase_edges = np.arange(-0.8,1.0,0.1)
    tavgs = OrderedDict()
    for i in range(len(phase_edges)-1):
        lo,hi = phase_edges[i],phase_edges[i+1]

        inds = np.nonzero(np.logical_and(phase_plus>lo,phase_plus<hi))[0]
        tavg = np.zeros(ds.nSamples)
        for i in inds:
            t=ds.read_trace(i)*1.0-ds.p_pretrig_mean[i]
            tavg+=t
        if len(inds)>30:
            tavg/=len(inds)
            tavgs[lo,hi]=tavg,len(inds)

    plt.figure()
    for (lo,hi),(tavg,npulses)  in tavgs.items():
        print("phase_plus={:.2f}, {:.2f}".format(np.mean([lo,hi]),np.dot(tavg[1:],ds.filter.filt_aterms[0])))

        plt.plot(tavg, label="{:.2f}, {:d}".format(np.mean([lo,hi]),npulses))
    plt.legend(title="phase_plus,npulses")
    plt.xlim(ds.nPresamples, ds.nPresamples+6)
    plt.ylim(-1,np.amax([tavg[ds.nPresamples+6] for tavg,npulses in tavgs.values()]))
    plt.grid(True)
    plt.xlabel("Sample Number")
    plt.ylabel("fb signal")
    plt.title("{}, chan {:d}\np_row_adjust = -(p_shift1+p_filt_phase)\nr=(p_rowcount+p_row_adjust*number_of_rows)\nr is subsample accurate arrival time estimate in row units".format(runstr,ds.channum))

def hist_ds(ds, edges, category={}):
    centers = 0.5*(edges[1:]+edges[:-1])
    counts, _ = np.histogram(ds.p_energy[ds.good(**category)],edges)
    return counts, centers

def hist_data(data, edges, category={}):
    centers = 0.5*(edges[1:]+edges[:-1])
    counts = np.zeros(len(centers))
    for ds in data:
        dscounts,_ = hist_ds(ds,edges,category=category)
        counts += dscounts
    return counts, centers
        
