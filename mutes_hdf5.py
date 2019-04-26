import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import numpy as np
import os
import h5py
import sys
import datetime as dt
import math
# -------- for ROOT file --------------
try:
    import ROOT
except ImportError:
    raise ValueError('ERROR: cannot import pyROOT')
ROOT.gROOT.SetBatch(1)

def get_name(runs,chs=[None]):
    name = "_run"+str(runs[0])
    if len(runs)>1: name+="_"+str(runs[-1])
    name += "_ch"+str(chs[0])
    if len(chs)>1: name+="_"+str(chs[-1])
    return name

def fill_1d(hist,hdf,chs,attr,cutlist):
    for ch,cut in zip(chs,cutlist):
        try:
            val=get_attr_float(hdf,ch,attr,cut)
            n=len(val)
            if n==0: continue
            w=np.ones(n,dtype=np.float64)
            hist.FillN(n,val,w)
        except:
            pass
            
def fill_2d(hist,hdf,chs,attr1,attr2,cutlist):
    for ch,cut in zip(chs,cutlist):
        try:
            val1=get_attr_float(hdf,ch,attr1,cut)
            val2=get_attr_float(hdf,ch,attr2,cut)
            n=len(val1)
            w=np.ones(n,dtype=np.float64)
            hist.FillN(n,val1,val2,w)
        except:
            pass

def mkroot_1d(hname,bins,hdf,chs,attr,cutlist):
    hist=ROOT.TH1F(hname,hname,bins[0],bins[1],bins[2])
    fill_1d(hist,hdf,chs,attr,cutlist)
    hist.Write()

def mkroot_2d(hname,bins,hdf,chs,attr1,attr2,cutlist):
    hist=ROOT.TH2F(hname,hname,bins[0],bins[1],bins[2],bins[3],bins[4],bins[5])
    fill_2d(hist,hdf,chs,attr1,attr2,cutlist)
    hist.Write()

def plot_hist(hdfs,chs,attr,edges,cutarray,ax=None):
    if ax==None:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)
    val,_=get_attr_sum(hdfs, chs, attr1=attr, attr2=None, cutarray=cutarray)
    H = ax.hist(val, bins=edges,histtype='step',alpha=0.7)
    return ax,H

def plot_hist2d(hdfs,chs,attr1,attr2,edges1,edges2,cutarray,ax=None):
    if ax==None:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)
    val1,val2=get_attr_sum(hdfs, chs, attr1, attr2, cutarray=cutarray)
    ax.hist2d(val1,val2, bins=[edges1,edges2],norm=LogNorm())
    return ax

def plot_scat(hdfs,chs,attr1,attr2,edges1,edges2,cutarray,ax=None):
    if ax==None:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)
    val1,val2=get_attr_sum(hdfs, chs, attr1, attr2, cutarray=cutarray)
    ax.scatter(val1,val2)
    return ax

#no need?
def plot_scat(val1,val2):
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(val1, val2)
    plt.show()
    return ax 

def get_filename(run,cut_pre=0,cut_post=0,DATADIR=None):
    if DATADIR==None:
        DATADIR = os.environ.get("MUTESDATADIR")
    if os.path.exists(DATADIR)==False:
        print "%s is missing"%DATADIR
        sys.exit()
    add='mass'
    if cut_pre > 0 or cut_post>0:
        add  += str("_pre%03d_post%03d" %(cut_pre,cut_post))
    filename= DATADIR + "/run" + str("%04d" %run) + "/run" + str("%04d" %run) + add + ".hdf5"
    if os.path.exists(filename)==False:
        print "%s is missing"%filename
        return ValueError
    return filename

def open_file(run,DATADIR=None):
    filename=get_filename(run,0,0,DATADIR)
    h5=h5py.File(filename,'r')
    return h5

def get_cut(h5,ch,cutlist):
    cutarray = True
    for cut in cutlist:
        tmpcut = get_cutarray(h5,ch,cut[0],cut[1],cut[2])
        cutarray= np.logical_and(cutarray,tmpcut)
    return cutarray

def get_cut_list(h5,chs,cutlist):
    alist= []
    for ch in chs:
        cutarray = True 
        for cut in cutlist:
            tmpcut = get_cutarray(h5,ch,cut[0],cut[1],cut[2])
            cutarray= np.logical_and(cutarray,tmpcut)
        alist.append(cutarray)
    return alist

def get_cutarray(h5,ch,attr,cut1,cut2):
    val=get_attr_np(h5,ch,attr)
    if cut2==None:
        return val == cut1
    cutarray1 = (val > cut1)
    cutarray2 = (val < cut2)
    return np.logical_and(cutarray1,cutarray2)

def get_attr_sum(hdfs, chs, attr1, attr2=None, cutarray=[]):
    val1=[]
    val2=[]
    for hdf in hdfs:
        for ch in chs:
            try:
                cut=get_cut(hdf,ch,cutarray)
                val1.extend(get_attr(hdf,ch,attr1,cut).tolist())
                if attr2==None: continue
                val2.extend(get_attr(hdf,ch,attr2,cut).tolist())
            except:
                pass
    return val1,val2

def get_attr(h5, ch, attr, cut=[None]):
    channel = 'chan' + str(ch)                        
    if channel not in h5: return None
    ds = h5[channel]
    if attr in ds:
        if cut[0]==None: return ds[attr]
        else: return ds[attr][cut]
    print "no attribute ", attr, " in ", h5.filename
    return None

def get_attr_np(h5, ch, attr, cut=[None]):
    return np.array(get_attr(h5,ch,attr,cut))

def get_attr_float(h5, ch, attr, cut=[None]):
    return np.array(get_attr(h5,ch,attr,cut),dtype=np.float64)


class mutes_hdf5():
    def __init__(self,runs=[],DATADIR=None):
        self.runs=runs
        self.hdf5_dict={}
        for run in self.runs:
            try:
                h5=open_file(run)
                self.hdf5_dict[run]=h5                
            except:
                continue
    def get_hdfs(self, runs):
        hdfs=[]
        for run in runs:
            if run in self.hdf5_dict:
                hdfs.append(self.hdf5_dict[run])
        return hdfs

    def get_attr(self,run, ch, attr):
        if run not in self.hdf5_dict:  return None
        h5=self.hdf5_dict[run]
        return get_attr(h5,ch,attr)

    def get_calib(self,run, ch, attr='p_filt_value_dc'):
        temp = self.get_attr(run,ch,'calibration')
        if temp==None: 
            print "no calibration for run", run, "ch", ch, attr
            return None, None
        if attr in temp:
            ph=temp[attr]['ph']
            name=temp[attr]['name']
            return ph.value.tolist(),name.value.tolist()
        return None, None

    def get_calib_peaks(self,run, ch, peaks=['CuKAlpha'], attr='p_filt_value_dc'):
        ph, name = self.get_calib(run,ch,attr)
        if ph==None: return None
        phs=[]
        for peak in peaks:
            try:
                index = name.index(peak)
            except ValueError:
                return ValueError        
            phs.append(ph[index])
        return phs

    def get_time_range(self,run, ch):
        temp = self.get_attr(run,ch,'timestamp')
        if temp==None:            return  [0,0]
        return [temp[0],temp[-1]]

    def get_time_err(self,run, ch):
        time = self.get_attr(run,ch,'timestamp')
        if time==None:            return  [0,0]
        mean = 0.5 * (time[0] + time[-1])
        err = 0.5 * (time[-1] - time[0])
        return dt.datetime.fromtimestamp(mean),dt.timedelta(seconds=err)

    def plot_calib_trend(self,runs,chs,peaks=['CuKAlpha'], attr='p_filt_value_dc'):
        pdfname="calib_trend"+get_name(runs,chs)+".pdf"
        print "printing...", pdfname
        with PdfPages(pdfname) as pdf:
            for ch in chs:
                time=[]
                err=[]
                phs=[]
                for run in runs:
                    t,e=self.get_time_err(run,ch)
                    ph=self.get_calib_peaks(run,ch,peaks,attr)
                    time.append(t)
                    err.append(e)
                    phs.append(ph)
                phs_t = np.array(phs).T.tolist()
                title="peak_position"+get_name(runs,[ch]) 
                plot_trend(time=time,etime=err,ylist=phs_t,names=peaks,title=title)
                pdf.savefig()
                plt.close()

    def plot_hist2d(self,runs,chs,attr1,attr2,edges1,edges2,cutarray,ax=None):
        hdfs=self.get_hdfs(runs)
        return plot_hist2d(hdfs,chs,attr1,attr2,edges1,edges2,cutarray,ax)

    def plot_hist(self,runs,chs,attr,edges,cutarray,ax=None):
        hdfs=self.get_hdfs(runs)
        return plot_hist(hdfs,chs,attr,edges,cutarray,ax)

    def plot_trend_diff(self,runs,chs,attr,cutarray):
        hdfs=self.get_hdfs(runs)
        for ch in chs:
            time=[]
            err=[]
            phs=[]
            figname="fig/trend_"+attr+get_name(runs,chs)+".png"
            plot_trend_diff(hdfs,[ch],attr,cutarray)
            print "printing...",figname
            plt.savefig(figname)
            plt.close()

    def plot_trend_multi(self,runs,chs,attrs,cutarray):
        hdfs=self.get_hdfs(runs)
        for ch in chs:
            time=[]
            err=[]
            phs=[]
            figname="fig/trend_run"+str(runs[0])+"_"+str(runs[-1])+"_ch"+str(ch)
            for attr in attrs:
                figname+="_"+attr
            figname+=".png"
            print "printing...",figname
            title="run"+str(runs[0])+"_"+str(runs[-1])+"_ch"+str(ch)
            plot_trend_multi(hdfs,[ch],attrs,cutarray,title)
            plt.savefig(figname)
            plt.close()

    def plot_average_pulse(self,runs,chs,attr='average_pulse',DIFF=False):
        pdfname="fig/"+attr.replace('/','_')+"_run"+str(runs[0])+"_"+str(runs[-1])
        if DIFF: pdfname+="_diff"
        pdfname+=".png"
        print "printing...", pdfname
        dx, dy=get_divisions(len(chs))
        fig = plt.figure()
        for index, ch in enumerate(chs):
            ax=fig.add_subplot(dx,dy,index+1)
            plt.title('chan'+str(ch))
            ref = self.get_attr(runs[0],ch,attr).value.tolist()
            if ref==None: continue
            index=range(0,len(ref))
            for run in runs:
                avg = self.get_attr(run,ch,attr).value.tolist()
                if avg==None: continue
                if DIFF: 
                    avg = [ y1 - y2 for y1, y2 in zip(avg,ref) ] 
                ax.plot(index,avg,label=str(run))
                ax.legend(fontsize=7)
        fig.tight_layout()
        plt.savefig(pdfname)
        plt.close()


def get_divisions(n):
    nS=math.sqrt(n)
    nc=int(math.ceil(nS))
    nr=int(nc-(nc*nc-n>nc-1))
    print n,nc,nr
    return nc, nr

def plot_trend(time,etime,ylist,names,title=''):
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    
    for y,name in zip(ylist,names):
        diffy=[ (tmp - y[0]) / y[0] * 100 for tmp in y ]
        ax1.errorbar(time,y,xerr=etime,fmt='.',label=name)
        ax2.errorbar(time,diffy,xerr=etime,fmt='.',label=name)
        
    ax1.set_title(title)
    ax1.set_ylabel('absolute')
    ax2.set_ylabel('relative[%]')        
    ax1.legend()
    #    ax.xticks(rotation=60)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    #    fig.autofmt_xdate()

def plot_trend_diff(hdfs,chs,attr,cutarray):
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    time,val=get_attr_sum(hdfs,chs,'timestamp',attr,cutarray)
    time=[dt.datetime.fromtimestamp(t) for t in time]
    mean=np.array(val).mean()
    diffy=[ (tmp - mean) for tmp in val ]

    ax1.scatter(time,val,marker='.')
    ax2.scatter(time,diffy,marker='.')
    ax1.set_xlim(min(time),max(time))
    ax2.set_xlim(min(time),max(time))
    ax1.set_ylabel('absolute')
    ax2.set_ylabel('relative')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%h%d'))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%h%d'))

    
def plot_trend_multi(hdfs,chs,attrs,cutarray,title=''):
    fig = plt.figure()
    n = len(attrs)
    for i,attr in enumerate(attrs):
        ax = fig.add_subplot(n,1,i+1)
        time,val=get_attr_sum(hdfs,chs,'timestamp',attr,cutarray)
        time=[dt.datetime.fromtimestamp(t) for t in time]
        mean=np.array(val).mean()
        ax.scatter(time,val,marker='.')
        ax.set_xlim(min(time),max(time))
        ax.set_ylabel(attr)
#        print max(time)-min(time), 
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%h%d\n%H:%M'))
        if i+1 != n: ax.set_xticklabels([])
        if i==0:     ax.set_title(title)
    fig.tight_layout()
    
