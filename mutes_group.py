""" mutes_ext.py is a function package for MUTES E62 group trigger analysis 

History: 
2018-08-16 ; ver 1.0; branched from mutes_util.py
2019-02-01 ; ver 1.1 a little bit modified, Hideyuki Tatsuno

"""

__version__ = '1.1'

import mass
import mutes_util as util
import numpy as np
import h5py
import sys


def calc_group_trigger_params(data,GRTINFO,forceNew=False,maxgrchs=4):
    # GRTINFO is group trigger channel map
    for ds in data:
        if not forceNew:
            print "ch%d group trigger ana skipped."%(ds.channum)
            try:
                print "loading secondary peak region data from hdf5 file."
                ds.sec_pr_sum      = ds.hdf5_group["sec_pr_sum"]
                ds.sec_pr_max      = ds.hdf5_group["sec_pr_max"]
                ds.sec_pr_mean     = ds.hdf5_group["sec_pr_mean"]
                ds.sec_pr_maxmean  = ds.hdf5_group["sec_pr_maxmean"]
                ds.sec_pr_meanmean = ds.hdf5_group["sec_pr_meanmean"]
                ds.sec_pr_maxmin   = ds.hdf5_group["sec_pr_maxmin"]
                ds.sec_pr_meanmin  = ds.hdf5_group["sec_pr_meanmin"]
                ds.sec_enemean     = ds.hdf5_group["sec_enemean"]
                ds.sec_enemax      = ds.hdf5_group["sec_enemax"]
                ds.sec_fdmean      = ds.hdf5_group["sec_fdmean"]
                ds.sec_fdmax       = ds.hdf5_group["sec_fdmax"]
            except:
                print "ch%d group trigger ana is not done yet."%(ds.channum)
                return False
            continue
        elif forceNew: print "ch%d group trig analysis ...."%(ds.channum)
            
        npl = ds.nPulses
        p_rows = np.array(ds.p_rowcount[:])
        grti = np.array(ds.p_grouptrig)# group trig ch info
        grtchs = np.zeros(maxgrchs,dtype=np.intc)# group hit channels
        gdict = util.get_grtdict(GRTINFO)
        grtchs = gdict.get(ds.channum)
        # primary flag
        prima = np.zeros(shape=grti.shape,dtype=bool)
        prim_ind = np.where(grti==-1)[0]# -1 is primary event, other number is a secondary one
        prima[prim_ind]=True
        p_prima_rows=p_rows[prim_ind]
        # initialize
        ds.sec_pr_sum      = np.zeros(npl,dtype=np.float32)
        ds.sec_pr_max      = np.zeros(npl,dtype=np.float32)
        ds.sec_pr_mean     = np.zeros(npl,dtype=np.float32)
        ds.sec_pr_maxmean  = np.zeros(npl,dtype=np.float32)
        ds.sec_pr_meanmean = np.zeros(npl,dtype=np.float32)
        ds.sec_pr_maxmin   = np.zeros(npl,dtype=np.float32)
        ds.sec_pr_meanmin  = np.zeros(npl,dtype=np.float32)
        # ds.sec_pr_meanmax  = np.zeros(npl,dtype=np.float32)      
        ds.sec_enemean     = np.zeros(npl,dtype=np.float32) # energy : mean
        ds.sec_enemax      = np.zeros(npl,dtype=np.float32) # energy : max
        ds.sec_fdmean      = np.zeros(npl,dtype=np.float32) # filtered value : mean
        ds.sec_fdmax       = np.zeros(npl,dtype=np.float32) # filtered value : max

        sene=[] # secondary energy
        sfil=[] # secondary filetered value        
        sprs=[]
        sprm=[]
        sprx=[]
        for j,ch in enumerate(grtchs):
            if ch==0 or data.channel.has_key(ch)==False: continue
            ds_sec = data.channel[ch]
            sec_grti = np.array(ds_sec.p_grouptrig)
            sec_index = np.where(sec_grti==ds.channum)[0]# secondary event for this primary ds
            # shold be len(sec_index) == len(prim_ind)
            # 20190201 Kokode kokeru?
            if not len(sec_index) == len(prim_ind):
                print "length of secondary index is strange.... ",len(prim_ind),len(sec_index)
                raise Exception
            sprs.append(np.array(ds_sec.p_peak_region_sum)[sec_index])
            sprx.append(np.array(ds_sec.p_peak_region_max)[sec_index])
            sprm.append(np.array(ds_sec.p_peak_region_mean)[sec_index])            
            sene.append(np.array(ds_sec.p_energy)[sec_index])
            sfil.append(np.array(ds_sec.p_filt_value)[sec_index])
        
        for i, p_ind in enumerate(prim_ind):
            ds.sec_pr_sum[p_ind]       = max([sprs[j][i] for j in xrange(len(sprs))])
            ds.sec_pr_max[p_ind]       = max([sprx[j][i] for j in xrange(len(sprx))])
            ds.sec_pr_mean[p_ind]      = max([sprm[j][i] for j in xrange(len(sprm))]) ## same as sec_pr_meanmax            
            ds.sec_pr_maxmean[p_ind]   = np.average([sprx[j][i] for j in xrange(len(sprx))])
            ds.sec_pr_meanmean[p_ind]  = np.average([sprm[j][i] for j in xrange(len(sprm))])
            ds.sec_pr_maxmin[p_ind]    = min([sprx[j][i] for j in xrange(len(sprx))])
            ds.sec_pr_meanmin[p_ind]   = min([sprm[j][i] for j in xrange(len(sprm))])
            # ds.sec_pr_meanmax[p_ind]   = max([sprm[j][i] for j in xrange(len(sprm))]) ## RH
            ds.sec_enemean[p_ind]      = np.average([sene[j][i] for j in xrange(len(sprx))])
            ds.sec_enemax[p_ind]       = np.amax([sene[j][i] for j in xrange(len(sprx))])
            ds.sec_fdmean[p_ind]       = np.average([sfil[j][i] for j in xrange(len(sprx))])
            ds.sec_fdmax[p_ind]        = np.amax([sfil[j][i] for j in xrange(len(sprx))])
            
        h5_sec_pr_sum     = ds.hdf5_group.require_dataset("sec_pr_sum",(npl,),dtype=np.float32)
        h5_sec_pr_sum[:]  = ds.sec_pr_sum
        h5_sec_pr_mean    = ds.hdf5_group.require_dataset("sec_pr_mean",(npl,),dtype=np.float32)
        h5_sec_pr_mean[:] = ds.sec_pr_mean
        h5_sec_pr_max     = ds.hdf5_group.require_dataset("sec_pr_max",(npl,),dtype=np.float32)
        h5_sec_pr_max[:]  = ds.sec_pr_max
        ### YI
        h5_sec_pr_maxmean = ds.hdf5_group.require_dataset("sec_pr_maxmean",(npl,),dtype=np.float32)
        h5_sec_pr_maxmean[:] = ds.sec_pr_maxmean
        h5_sec_pr_meanmean = ds.hdf5_group.require_dataset("sec_pr_meanmean",(npl,),dtype=np.float32)
        h5_sec_pr_meanmean[:] = ds.sec_pr_meanmean
        h5_sec_pr_maxmin = ds.hdf5_group.require_dataset("sec_pr_maxmin",(npl,),dtype=np.float32)
        h5_sec_pr_maxmin[:] = ds.sec_pr_maxmin
        h5_sec_pr_meanmin = ds.hdf5_group.require_dataset("sec_pr_meanmin",(npl,),dtype=np.float32)
        h5_sec_pr_meanmin[:] = ds.sec_pr_meanmin
        ## RH
        # h5_sec_pr_meanmax = ds.hdf5_group.require_dataset("sec_pr_meanmax",(npl,),dtype=np.float32)
        # h5_sec_pr_meanmax[:] = ds.sec_pr_meanmax
        ## _RH
        ### _YI
        ### SY
        h5_sec_enemean = ds.hdf5_group.require_dataset("sec_enemean",(npl,),dtype=np.float32)
        h5_sec_enemean[:] = ds.sec_enemean
        h5_sec_enemax = ds.hdf5_group.require_dataset("sec_enemax",(npl,),dtype=np.float32)
        h5_sec_enemax[:] = ds.sec_enemax
        h5_sec_fdmean = ds.hdf5_group.require_dataset("sec_fdmean",(npl,),dtype=np.float32)
        h5_sec_fdmean[:] = ds.sec_fdmean
        h5_sec_fdmax = ds.hdf5_group.require_dataset("sec_fdmax",(npl,),dtype=np.float32)
        h5_sec_fdmax[:] = ds.sec_fdmax
        ### _SY
        
        ds.sec_pr_sum       = ds.hdf5_group["sec_pr_sum"]
        ds.sec_pr_max       = ds.hdf5_group["sec_pr_max"]
        ds.sec_pr_mean      = ds.hdf5_group["sec_pr_mean"]
        ds.sec_pr_maxmean   = ds.hdf5_group["sec_pr_maxmean"]
        ds.sec_pr_meanmean  = ds.hdf5_group["sec_pr_meanmean"]
        ds.sec_pr_maxmin    = ds.hdf5_group["sec_pr_maxmin"]
        ds.sec_pr_meanmin   = ds.hdf5_group["sec_pr_meanmin"]
        # ds.sec_pr_meanmax  = ds.hdf5_group["sec_pr_meanmax"]
        ds.sec_enemean      = ds.hdf5_group["sec_enemean"]
        ds.sec_enemax       = ds.hdf5_group["sec_enemax"]
        ds.sec_fdmean       = ds.hdf5_group["sec_fdmean"]
        ds.sec_fdmax        = ds.hdf5_group["sec_fdmax"]
        
    return True

def get_grouptrig_pulse_info(data, cutarray, cut_category,GRTINFO,getMax=100,chanlist=[129]):
    # group trigger channel map
    # GRTINFO = "../csv/grptrig_twocol.txt"
    GRTINFO = GRTINFO
    maxgrchs=4
    pulse_records_array_all = []
    peak_region_max_all = []
    peak_region_mean_all = []
    peak_region_sum_all = []
    energy_all = []
    pretrig_mean_all = []
    pretrig_rms_all = []
    peak_value_all = []

    print "chanlist = ", chanlist

    for ds in data:
#        for c, chan in enumerate(chanlist[0]):
        for k, chan in enumerate(chanlist):
            pulse_records_array_onechan = []
            peak_region_max_onechan = []
            peak_region_mean_onechan = []
            peak_region_sum_onechan = []
            energy_onechan = []
            pretrig_mean_onechan = []
            pretrig_rms_onechan = []
            peak_value_onechan = []
            if not ds.channum == chan:
                continue
            # if len(cutarray[c]) > 0:
            #     energy_array = np.array(ds.p_energy)
            #     sec_pr_mean_array = np.array(ds.sec_pr_mean) 
            #     if cut_category == "sec_pr_mean":
            #         cut_sec_pr_mean_array = sec_pr_mean_array[cutarray[c]]
            if ds.good().sum() == 0:
                print "Chan%d has no good pulses. Select (an)other channel(s)"%(ds.channum)
                continue

            npl = ds.nPulses
            grti = np.array(ds.p_grouptrig)# group trig ch info
            grtchs = np.zeros(maxgrchs,dtype=np.intc)# group hit channels
            gdict = util.get_grtdict(GRTINFO)
            grtchs = gdict.get(ds.channum)
            # primary flag
            p_rows = np.array(ds.p_rowcount[:])
            prima = np.zeros(shape=grti.shape,dtype=bool)
            #prim_ind = np.where(grti==-1)[0]
            prim_ind = np.where(grti==-1)[0]
            prima[prim_ind]=True
            p_prima_rows=p_rows[prim_ind]
            # initialize
            prima_cut = np.logical_and(prima,cutarray[k])
            prima_cut_ind = np.where(prima_cut==True)[0]
            prima_cut2 = np.logical_and(prima_cut, ds.good())
            prima_cut_ind2 = np.where(prima_cut2==True)[0]
            print "start getting pulses info / primary chan = ", chan
            print "total primary pulses =", len(prim_ind)
            print "total cut dataset =", len(prima_cut_ind2)
            if not len(grtchs) > 5:
                print "Num. of plot pulses =", getMax
                #print "Num. of plot dataset = ", int(getMax / (len(grtchs)+1))
            ev = 0
            ev2 = 0
            for i, oneidx in enumerate(prim_ind):

                ## for QL
                # if len(grtchs) > 5 and len(prima_cut_ind) > 5:
                #     if i > 3: print "break"; break
                ### 
                if not oneidx in prima_cut_ind2: 
                    #print "continue oneidx=",oneidx
                    continue
                # if ds.good()[oneidx] == False: continue
                # # if len(cutarray[c]) > 0:
                # #     if cut_category == "sec_pr_mean":
                # #         if not sec_pr_mean_array[oneidx] in cut_sec_pr_mean_array: continue
                ev2 += 1
                pulse_records_dict = {}
                peak_region_max_dict = {}
                peak_region_mean_dict = {}
                peak_region_sum_dict = {}
                energy_dict = {}
                pretrig_mean_dict = {}
                pretrig_rms_dict = {}
                peak_value_dict = {}

                pulse_records = []
                
                peak_region_max_dict[chan] = ds.p_peak_region_max[oneidx]
                peak_region_mean_dict[chan] = ds.p_peak_region_mean[oneidx]
                peak_region_sum_dict[chan] = ds.p_peak_region_sum[oneidx]
                energy_dict[chan] = ds.p_energy[oneidx]
                pretrig_mean_dict[chan] = ds.p_pretrig_mean[oneidx]
                pretrig_rms_dict[chan] = ds.p_pretrig_rms[oneidx]
                peak_value_dict[chan] = ds.p_peak_value[oneidx]

                for ii in xrange(ds.nSamples):
                    pulse_records.append(ds.pulse_records.datafile[oneidx][ii]-ds.p_pretrig_mean[oneidx])
                pulse_records = np.array(pulse_records)
                pulse_records_dict[chan] = pulse_records
                ev += 1
                for j,ch in enumerate(grtchs):

                    if data.channel.has_key(ch)==False:
                        #print "ch%d does not exist"
                        continue
                    ev += 1
                    ds_sec = data.channel[ch]
                    sec_grti = np.array(ds_sec.p_grouptrig)
                    seconda = np.zeros(shape=sec_grti.shape,dtype=bool)
                    sec_index = np.where(sec_grti==ds.channum)[0] # secondary for this ds
                    # seconda[sec_index]=True
                    # # initialize
                    # seconda_cut = np.logical_and(seconda,cutarray[c])
                    # seconda_cut_ind = np.where(seconda_cut==True)[0] 
                    # shold be len(sec_index) == len(prim_ind)
                    if not len(sec_index) == len(prim_ind):
                        print "length of secondary index is strange.... ",len(prim_ind),len(sec_index)
                        sys.exit()
                    # if not len(seconda_cut_ind) == len(prima_cut_ind):
                    #     print "length of secondary index is strange.... "
                    #     sys.exit()
                        # continue
                    peak_region_max_dict[ch] = ds_sec.p_peak_region_max[sec_index[i]]
                    peak_region_mean_dict[ch] = ds_sec.p_peak_region_mean[sec_index[i]]
                    peak_region_sum_dict[ch] = ds_sec.p_peak_region_sum[sec_index[i]]
                    energy_dict[ch] = ds_sec.p_energy[sec_index[i]]
                    pretrig_mean_dict[ch] = ds_sec.p_pretrig_mean[sec_index[i]]
                    pretrig_rms_dict[ch] = ds_sec.p_pretrig_rms[sec_index[i]]
                    peak_value_dict[ch] = ds_sec.p_peak_value[sec_index[i]]

                    pulse_records_array_sec = [] # to store secondary pulse records array

                    for jj in xrange(ds_sec.nSamples):
                        pulse_records_array_sec.append(ds_sec.pulse_records.datafile[sec_index[i]][jj]-ds_sec.p_pretrig_mean[sec_index[i]])
                    pulse_records_array_sec = np.array(pulse_records_array_sec)
                    pulse_records_dict[ch] = pulse_records_array_sec

                    # peak_region_max_dict[ch] = ds_sec.p_peak_region_max[seconda_cut_ind[i]]
                    # peak_region_mean_dict[ch] = ds_sec.p_peak_region_mean[seconda_cut_ind[i]]
                    # peak_region_sum_dict[ch] = ds_sec.p_peak_region_sum[seconda_cut_ind[i]]
                    # energy_dict[ch] = ds_sec.p_energy[seconda_cut_ind[i]]
                    # pretrig_mean_dict[ch] = ds_sec.p_pretrig_mean[seconda_cut_ind[i]]
                    # pretrig_rms_dict[ch] = ds_sec.p_pretrig_rms[seconda_cut_ind[i]]

                    # pulse_records_array_sec = [] # to store secondary pulse records array

                    # for jj in xrange(ds_sec.nSamples):
                    #     pulse_records_array_sec.append(ds_sec.pulse_records.datafile[seconda_cut_ind[i]][jj]-ds_sec.p_pretrig_mean[seconda_cut_ind[i]])
                    # pulse_records_array_sec = np.array(pulse_records_array_sec)
                    # pulse_records_dict[ch] = pulse_records_array_sec


                pulse_records_array_onechan.append(pulse_records_dict)
                peak_region_max_onechan.append(peak_region_max_dict)
                peak_region_mean_onechan.append(peak_region_mean_dict)
                peak_region_sum_onechan.append(peak_region_sum_dict)
                energy_onechan.append(energy_dict)
                pretrig_mean_onechan.append(pretrig_mean_dict)
                pretrig_rms_onechan.append(pretrig_rms_dict)
                peak_value_onechan.append(peak_value_dict)
                if not len(grtchs) > 5: print "get pulse info %d/%d"%(ev2,int(getMax / len(pulse_records_dict.keys())))
                else: print "get pulse info %d/%d"%(ev2,len(prima_cut_ind2))
                if ev == getMax: print "chan%d done"%(chan); break

            pulse_records_array_all.append(pulse_records_array_onechan)
            peak_region_max_all.append(peak_region_max_onechan)
            peak_region_mean_all.append(peak_region_mean_onechan)
            peak_region_sum_all.append(peak_region_sum_onechan)
            energy_all.append(energy_onechan)
            pretrig_mean_all.append(pretrig_mean_onechan)
            pretrig_rms_all.append(pretrig_rms_onechan)
            peak_value_all.append(peak_value_onechan)

    return pulse_records_array_all, peak_region_max_all, peak_region_mean_all, peak_region_sum_all, energy_all, pretrig_mean_all, pretrig_rms_all, peak_value_all


