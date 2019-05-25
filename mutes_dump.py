""" mutes_dump.py is a function package to dumping into a ROOT file for MUTES E62 

History: 
2018-08-16 ; ver 1.0; branched from mutes_ana.py, T.H.
2019-02-19 ; ver 1.1; HT
2019-02-22 ; ver 1.2; HT modified names of numpy arrays and Fill command
"""

__version__ = '1.2'


import time
import mass
import mutes_ext as ext
import mutes_util as util
import numpy as np
import os

from ROOT import gROOT
gROOT.SetBatch(1)
import ROOT

ROOTDIR="%s/dumproot"%(os.environ['MUTESDATADIR'])


def root_hist1d(values,inds=[],name="h1",nbin=10,minx=0,maxx=10,title="",xtitle="",ytitle=""):
    h = ROOT.TH1F(name,name,nbin,minx,maxx)
    if len(inds)==0: inds=np.arange(len(values))
    n = len(inds)
    h.FillN(n,np.array(values,dtype=np.float64)[inds],np.ones(n,dtype=np.float64))
    h.Sumw2()
    h.SetTitle(title)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    return h


def dump_ROOT(data,fout="hoge.root",EXTTRIG=True, GRTRIG=True, dumppulse=False):
    
    start = time.time()

    # --- TFile ---
    f = ROOT.TFile(fout,"recreate")
    cdhist=f.mkdir("hist")
    f.cd()
    print "...dump_ROOT start, save to %s"%fout
    # --- TTree ---
    cname, ctitle = 'ct', 'common tree'
    ct = ROOT.TTree(cname,ctitle)
    ctfill = ct.Fill
    ctbranch = ct.Branch
    bc_run              = np.zeros(1,dtype=np.intc)
    bc_ch               = np.zeros(1,dtype=np.intc)
    bc_col              = np.zeros(1,dtype=np.intc)
    bc_row              = np.zeros(1,dtype=np.intc)
    bc_ncols            = np.zeros(1,dtype=np.intc)
    bc_nrows            = np.zeros(1,dtype=np.intc)
    bc_npulses          = np.zeros(1,dtype=np.intc)
    bc_nsamples         = np.zeros(1,dtype=np.intc)
    bc_npresamples      = np.zeros(1,dtype=np.intc)
    bc_row_timebase     = np.zeros(1,dtype=np.float64)
    bc_timebase         = np.zeros(1,dtype=np.float64)
    bc_timestamp_offset = np.zeros(1,dtype=np.float64)
    ctbranch('run',     bc_run     ,'run/I')
    ctbranch('ch',      bc_ch      ,'ch/I')
    ctbranch('col',     bc_col     ,'col/I')
    ctbranch('row',     bc_row     ,'row/I')
    ctbranch('ncols',   bc_ncols   ,'ncols/I')
    ctbranch('nrows',   bc_nrows   ,'nrows/I')
    # common parameters
    ds0 = data.first_good_dataset
    run              = ds0.runnum
    timebase         = ds0.timebase
    nsamples         = ds0.nSamples
    npresamples      = ds0.nPresamples
    timestamp_offset = ds0.timestamp_offset
    row_timebase     = timebase/float(ds0.number_of_rows)
    number_of_rows   = ds0.number_of_rows
    frame_time       = row_timebase*number_of_rows

    # --- TTree ---
    pname, ptitle = 'chanall', 'pulse tree'
    pt = ROOT.TTree(pname,ptitle)
    ptfill = pt.Fill
    ptbranch = pt.Branch
    bp_run             = np.zeros(1,dtype=np.intc)
    bp_ch              = np.zeros(1,dtype=np.intc)
    bp_ev              = np.zeros(1,dtype=np.intc)
    bp_good            = np.zeros(1,dtype=bool)
    bp_prim            = np.zeros(1,dtype=bool)
    bp_jbrsc           = np.zeros(1,dtype=bool)
    bp_grt             = np.zeros(1,dtype=np.intc)
    bp_shift1          = np.zeros(1,dtype=bool)
    bp_filt_phase      = np.zeros(1,dtype=np.float32)
    bp_filt_value      = np.zeros(1,dtype=np.float32)
    bp_filt_value_dc   = np.zeros(1,dtype=np.float32)
    bp_filt_value_phc  = np.zeros(1,dtype=np.float32)
    bp_filt_value_tdc  = np.zeros(1,dtype=np.float32)
    bp_energy          = np.zeros(1,dtype=np.float32)
    bp_energy_dc       = np.zeros(1,dtype=np.float32)
    bp_energy_phc      = np.zeros(1,dtype=np.float32)
    bp_min_value       = np.zeros(1,dtype=np.float32)
    bp_peak_index      = np.zeros(1,dtype=np.intc)
    bp_peak_time       = np.zeros(1,dtype=np.float32)
    bp_peak_value      = np.zeros(1,dtype=np.float32)
    bp_postpeak_deriv  = np.zeros(1,dtype=np.float32)
    bp_pretrig_mean    = np.zeros(1,dtype=np.float32)
    bp_pretrig_rms     = np.zeros(1,dtype=np.float32)
    bp_promptness      = np.zeros(1,dtype=np.float32)
    bp_pulse_average   = np.zeros(1,dtype=np.float32)
    bp_pulse_rms       = np.zeros(1,dtype=np.float32)
    bp_rise_time       = np.zeros(1,dtype=np.float32)
    bp_timestamp       = np.zeros(1,dtype=np.float64)
    bp_rowcount        = np.zeros(1,dtype=np.float64)
    bp_rowp            = np.zeros(1,dtype=np.float64)
    bp_rown            = np.zeros(1,dtype=np.float64)
    # --- cut paramaters ---
    bp_peak_region_max  = np.zeros(1,dtype=np.float32)                        
    bp_peak_region_mean = np.zeros(1,dtype=np.float32)
    bp_pre_region_sum   = np.zeros(1,dtype=np.float32)
    bp_jbr_region_sum   = np.zeros(1,dtype=np.float32)
    if EXTTRIG:
        bp_beamOn       = np.zeros(1,dtype=bool)
        bp_rows_after_last_external_trigger_nrp = np.zeros(1,dtype=np.int64)
        bp_rows_until_next_external_trigger_nrp = np.zeros(1,dtype=np.int64)
        bp_rows_after_last_external_trigger_nrn = np.zeros(1,dtype=np.int64)
        bp_rows_until_next_external_trigger_nrn = np.zeros(1,dtype=np.int64)
    if GRTRIG:
        bp_sprmc           = np.zeros(1,dtype=bool)
        bp_sec_pr_maxmean  = np.zeros(1,dtype=np.float32)            
        bp_sec_pr_meanmean = np.zeros(1,dtype=np.float32)            
        bp_sec_pr_mean     = np.zeros(1,dtype=np.float32)            
        bp_sec_enemean     = np.zeros(1,dtype=np.float32)            
        bp_sec_enemax      = np.zeros(1,dtype=np.float32)            

    # --- tree branch address
    ptbranch('run',            bp_run,             'run/I')
    ptbranch('ch',             bp_ch,              'ch/I')
    ptbranch('ev',             bp_ev,              'ev/I')
    ptbranch('good',           bp_good,            'good/O')
    ptbranch('primary',        bp_prim,            'primary/O')
    ptbranch('jbrsc',          bp_jbrsc,           'jbrsc/O')
    ptbranch('grouptrigger',   bp_grt,             'grouptrigger/I')
    ptbranch('shift1',         bp_shift1,          'shift1/O')
    ptbranch('filt_phase',     bp_filt_phase,      'filt_phase/F')
    ptbranch('filt_value',     bp_filt_value,      'filt_value/F')
    ptbranch('filt_value_dc',  bp_filt_value_dc,   'filt_value_dc/F')
    ptbranch('filt_value_phc', bp_filt_value_phc,  'filt_value_phc/F')
    ptbranch('filt_value_tdc', bp_filt_value_tdc,  'filt_value_tdc/F')
    ptbranch('energy',         bp_energy,          'energy/F')
    ptbranch('energy_dc',      bp_energy_dc,       'energy_dc/F')
    ptbranch('energy_phc',     bp_energy_phc,      'energy_phc/F')
    ptbranch('min_value',      bp_min_value,       'min_value/F')
    ptbranch('peak_index',     bp_peak_index,      'peak_index/I')
    ptbranch('peak_time',      bp_peak_time,       'peak_time/F')
    ptbranch('peak_value',     bp_peak_value,      'peak_value/F')
    ptbranch('postpeak_deriv', bp_postpeak_deriv,  'postpeak_deriv/F')
    ptbranch('pretrig_mean',   bp_pretrig_mean,    'pretrig_mean/F')
    ptbranch('pretrig_rms',    bp_pretrig_rms,     'pretrig_rms/F')
    ptbranch('promptness',     bp_promptness,      'promptness/F')
    ptbranch('pulse_average',  bp_pulse_average,   'pulse_average/F')
    ptbranch('pulse_rms',      bp_pulse_rms,       'pulse_rms/F')
    ptbranch('rise_time',      bp_rise_time,       'rise_time/F')
    ptbranch('timestamp',      bp_timestamp,       'timestamp/D')
    ptbranch('rowcount',       bp_rowcount,        'rowcount/D')
    ptbranch('rowp',           bp_rowp,            'rowp/D')
    ptbranch('rown',           bp_rown,            'rown/D')
    # --- cut paramaters 
    ptbranch('peak_region_max',  bp_peak_region_max,  'peak_region_max/F')
    ptbranch('peak_region_mean', bp_peak_region_mean, 'peak_region_mean/F')
    ptbranch('pre_region_sum',   bp_pre_region_sum,   'pre_region_sum/F')
    ptbranch('jbr_region_sum',   bp_jbr_region_sum,   'jbr_region_sum/F')
    if EXTTRIG:
        ptbranch('beam',                bp_beamOn,                               'beam/O')
        ptbranch('row_after_extrig_nrp',bp_rows_after_last_external_trigger_nrp, 'row_after_extrig_nrp/L')
        ptbranch('row_next_extrig_nrp', bp_rows_until_next_external_trigger_nrp, 'row_next_extrig_nrp/L')
        ptbranch('row_after_extrig_nrn',bp_rows_after_last_external_trigger_nrn, 'row_after_extrig_nrn/L')
        ptbranch('row_next_extrig_nrn', bp_rows_until_next_external_trigger_nrn, 'row_next_extrig_nrn/L')
    if GRTRIG:
        ptbranch('sprmc',           bp_sprmc,           'sprmc/O')
        ptbranch('sec_pr_maxmean',  bp_sec_pr_maxmean,  'sec_pr_maxmean/F')
        ptbranch('sec_pr_meanmean', bp_sec_pr_meanmean, 'sec_pr_meanmean/F')
        ptbranch('sec_pr_mean',     bp_sec_pr_mean,     'sec_pr_mean/F')
        ptbranch('sec_enemean',     bp_sec_enemean,     'sec_enemean/F')
        ptbranch('sec_enemax',      bp_sec_enemax,      'sec_enemax/F')

    # loop start
    for ds in data:
        dschan=ds.channum
        print "channel %d start.... for %.3f (sec)"%(dschan, (time.time() - start))       
        # ---------------------------------------------
        bc_run[0]              = run
        bc_ch[0]               = dschan
        bc_col[0]              = ds.column_number
        bc_row[0]              = ds.row_number
        bc_ncols[0]            = ds.number_of_columns
        bc_nrows[0]            = ds.number_of_rows
        bc_npulses[0]          = ds.nPulses
        bc_nsamples[0]         = ds.nSamples
        bc_npresamples[0]      = ds.nPresamples
        bc_row_timebase[0]     = ds.timebase/float(ds.number_of_rows)
        bc_timebase[0]         = ds.timebase
        bc_timestamp_offset[0] = ds.timestamp_offset
        ctfill()
        # ---------------------------------------------            
        # --- numpy array objects for fast calculation
        np_good          =np.array(ds.p_goodflag)
        np_prim          =np.array(ds.p_prime)
        np_jbrsc         =np.array(ds.p_jbrsc)
        np_grt           =np.array(ds.p_grouptrig)
        np_shift1        =np.array(ds.p_shift1)
        np_filt_phase    =np.array(ds.p_filt_phase)
        np_filt_value    =np.array(ds.p_filt_value)
        np_filt_value_dc =np.array(ds.p_filt_value_dc)
        np_filt_value_phc=np.array(ds.p_filt_value_phc)
        np_filt_value_tdc=np.array(ds.p_filt_value_tdc)
        np_energy        =np.array(ds.p_energy)
        np_energy_dc     =np.zeros(shape=np_filt_value_dc.shape,dtype=np.float32)
        np_energy_phc    =np.zeros(shape=np_filt_value_phc.shape,dtype=np.float32)
        attr="p_filt_value_dc"
        if ds.calibration.has_key(attr):
            np_energy_dc = ds.calibration[attr].ph2energy(np_filt_value_dc)
        else:
            print "no calibration %s in ch%d"%(attr,dschan)
        attr="p_filt_value_phc"
        if ds.calibration.has_key(attr):
            np_energy_phc = ds.calibration[attr].ph2energy(np_filt_value_phc)
        else:
            print "no calibration %s in ch%d"%(attr,dschan)
        np_min_value     =np.array(ds.p_min_value)
        np_peak_index    =np.array(ds.p_peak_index)
        np_peak_time     =np.array(ds.p_peak_time)
        np_peak_value    =np.array(ds.p_peak_value)
        np_postpeak_deriv=np.array(ds.p_postpeak_deriv)
        np_pretrig_mean  =np.array(ds.p_pretrig_mean)
        np_pretrig_rms   =np.array(ds.p_pretrig_rms)
        np_promptness    =np.array(ds.p_promptness)
        np_pulse_average =np.array(ds.p_pulse_average)
        np_pulse_rms     =np.array(ds.p_pulse_rms)
        np_rise_time     =np.array(ds.p_rise_time)
        np_timestamp     =np.array(ds.p_timestamp)
        np_rowcount      =np.array(ds.p_rowcount)
        np_rown          =np.array(ds.p_rown)
        np_rowp          =np.array(ds.p_rowp)
        # --- cut paramaters 
        np_peak_region_max =np.array(ds.p_peak_region_max)
        np_peak_region_mean=np.array(ds.p_peak_region_mean)
        np_pre_region_sum  =np.array(ds.p_pre_region_sum)
        np_jbr_region_sum  =np.array(ds.p_jbr_region_sum)
        # --- group trigger ----
        if GRTRIG:
            np_sprmc          =np.array(ds.p_sprmc)
            np_sec_pr_maxmean =np.array(ds.sec_pr_maxmean)
            np_sec_pr_meanmean=np.array(ds.sec_pr_meanmean)
            np_sec_pr_mean    =np.array(ds.sec_pr_mean)
            np_sec_enemean    =np.array(ds.sec_enemean)
            np_sec_enemax     =np.array(ds.sec_enemax)
        # --- external trigger ----
        if EXTTRIG:
            np_beamOn                              =np.array(ds.p_beamflag)
            np_rows_after_last_external_trigger_nrp=np.array(ds.rows_after_last_external_trigger_nrp)
            np_rows_until_next_external_trigger_nrp=np.array(ds.rows_until_next_external_trigger_nrp)
            np_rows_after_last_external_trigger_nrn=np.array(ds.rows_after_last_external_trigger_nrn)
            np_rows_until_next_external_trigger_nrn=np.array(ds.rows_until_next_external_trigger_nrn)
        
        cdhist.cd()
        nbin = 25000
        minx = 0.
        maxx = 25000.
        hout_name = "%s_gpj%d_ch%d"%(util.hpht_phc,run,dschan)
        gp  = np.logical_and(np_good,np_prim)
        gpj = np.logical_and(gp,np_jbrsc)
        gpj_ind=np.where(gpj==True)[0]
        hout = root_hist1d(np_filt_value_phc,gpj_ind,hout_name,nbin,minx,maxx,
                           title="run%04d ch%d filt_value_phc gpj"%(run,dschan),
                           xtitle="filt_value_phc [ch]",ytitle="Counts / ch")
        hout.Write()
        hout.Delete()
        if EXTTRIG and GRTRIG:
            hout_onsprmon_name = "%s_onsprmon%d_ch%d"%(util.hpht_phc,run,dschan)#  beamon && sprmcon
            onsprmon = np.logical_and(np_beamOn,np_sprmc)
            gpj_onsprmon = np.logical_and(gpj,onsprmon)
            gpj_onsprmon_ind = np.where(gpj_onsprmon==True)[0]
            hout_onsprmon = root_hist1d(np_filt_value_phc,gpj_onsprmon_ind,hout_onsprmon_name,nbin,minx,maxx,
                                        title="run%04d ch%d filt_value_phc gpj_onsprmon"%(run,dschan),
                                        xtitle="filt_value_phc [ch]",ytitle="Counts / ch")
            hout_onsprmon.Write()
            hout_onsprmon.Delete()
            hout_offsprmon_name = "%s_offsprmon%d_ch%d"%(util.hpht_phc,run,dschan)#  beamon && sprmcon
            offsprmon = np.logical_and(~np_beamOn,np_sprmc)
            gpj_offsprmon = np.logical_and(gpj,offsprmon)
            gpj_offsprmon_ind = np.where(gpj_offsprmon==True)[0]
            hout_offsprmon = root_hist1d(np_filt_value_phc,gpj_offsprmon_ind,hout_offsprmon_name,nbin,minx,maxx,
                                        title="run%04d ch%d filt_value_phc gpj_offsprmon"%(run,dschan),
                                        xtitle="filt_value_phc [ch]",ytitle="Counts / ch")
            hout_offsprmon.Write()
            hout_offsprmon.Delete()
        f.cd()    

        # tree fill
        for i in xrange(ds.nPulses):
            bp_run[0]             = run
            bp_ch[0]              = dschan
            bp_ev[0]              = i
            bp_good[0]            = np_good[i]          
            bp_grt[0]             = np_grt[i]                     
            bp_prim[0]            = np_prim[i]
            bp_jbrsc[0]           = np_jbrsc[i]
            bp_shift1[0]          = np_shift1[i]                  
            bp_filt_phase[0]      = np_filt_phase[i]              
            bp_filt_value[0]      = np_filt_value[i]              
            bp_filt_value_dc[0]   = np_filt_value_dc[i]          
            bp_filt_value_phc[0]  = np_filt_value_phc[i]          
            bp_filt_value_tdc[0]  = np_filt_value_tdc[i]          
            bp_energy[0]          = np_energy[i]                  
            bp_energy_dc[0]       = np_energy_dc[i]               
            bp_energy_phc[0]      = np_energy_phc[i]              
            bp_min_value[0]       = np_min_value[i]               
            bp_peak_index[0]      = np_peak_index[i]              
            bp_peak_time[0]       = np_peak_time[i]               
            bp_peak_value[0]      = np_peak_value[i]              
            bp_postpeak_deriv[0]  = np_postpeak_deriv[i]          
            bp_pretrig_mean[0]    = np_pretrig_mean[i]            
            bp_pretrig_rms[0]     = np_pretrig_rms[i]             
            bp_promptness[0]      = np_promptness[i]              
            bp_pulse_average[0]   = np_pulse_average[i]           
            bp_pulse_rms[0]       = np_pulse_rms[i]               
            bp_rise_time[0]       = np_rise_time[i]               
            bp_timestamp[0]       = np_timestamp[i]               
            bp_rowcount[0]        = np_rowcount[i]                
            bp_rown[0]            = np_rown[i]                 
            bp_rowp[0]            = np_rowp[i]                 
            bp_peak_region_max[0] = np_peak_region_max[i]
            bp_peak_region_mean[0]= np_peak_region_mean[i]
            bp_pre_region_sum[0]  = np_pre_region_sum[i]
            bp_jbr_region_sum[0]  = np_jbr_region_sum[i]
            if EXTTRIG:
                bp_beamOn[0]      = np_beamOn[i]
                bp_rows_after_last_external_trigger_nrp[0] = np_rows_after_last_external_trigger_nrp[i]
                bp_rows_until_next_external_trigger_nrp[0] = np_rows_until_next_external_trigger_nrp[i]
                bp_rows_after_last_external_trigger_nrn[0] = np_rows_after_last_external_trigger_nrn[i]
                bp_rows_until_next_external_trigger_nrn[0] = np_rows_until_next_external_trigger_nrn[i]
            if GRTRIG:
                bp_sprmc[0]           = np_sprmc[i]                    
                bp_sec_pr_maxmean[0]  = np_sec_pr_maxmean[i]
                bp_sec_pr_meanmean[0] = np_sec_pr_meanmean[i]
                bp_sec_pr_mean[0]     = np_sec_pr_mean[i]
                bp_sec_enemean[0]     = np_sec_enemean[i]
                bp_sec_enemax[0]      = np_sec_enemax[i]
            ptfill()

            
    pt.Write()
    ct.Write()

    f.Close()
    print "End of dump to ROOT for %.3f sec"%((time.time() - start))
    print "-------------------------------------------------"
