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

#ROOTDIR="/Volumes/MUTES_HD/root"
ROOTDIR="%s/dumproot"%(os.environ['MUTESDATADIR'])

def dump_ROOT_2019(data,fout,EXTTRIG=True, GRTRIG=True, dumppulse=False):
    from ROOT import gROOT
    gROOT.SetBatch(1)
    import ROOT

    start = time.time()

    # --- TFile ---
    f = ROOT.TFile(fout,"recreate")
    print "...dump_ROOT start, save to %s"%fout
    # --- TTree ---
    cname, ctitle = 'ct', 'common tree'
    ct = ROOT.TTree(cname,ctitle)
    ctfill = ct.Fill
    ctbranch = ct.Branch
    bc_ch                 = np.zeros(1,dtype=np.intc)
    bc_col                = np.zeros(1,dtype=np.intc)
    bc_row                = np.zeros(1,dtype=np.intc)
    bc_ncols              = np.zeros(1,dtype=np.intc)
    bc_nrows              = np.zeros(1,dtype=np.intc)
    bc_npulses            = np.zeros(1,dtype=np.intc)
    bc_nsamples           = np.zeros(1,dtype=np.intc)
    bc_npresamples        = np.zeros(1,dtype=np.intc)
    bc_row_timebase       = np.zeros(1,dtype=np.float64)
    bc_timebase           = np.zeros(1,dtype=np.float64)
    bc_timestamp_offset   = np.zeros(1,dtype=np.float64)
    ctbranch('ch',      bc_ch      ,'channel number/I')
    ctbranch('col',     bc_col     ,'column number/I')
    ctbranch('row',     bc_row     ,'row number/I')
    ctbranch('ncols',   bc_ncols   ,'n of columns/I')
    # common parameters
    ds0 = data.first_good_dataset
    timebase    = ds0.timebase
    nsamples    = ds0.nSamples
    npresamples = ds0.nPresamples
    timestamp_offset = ds0.timestamp_offset
    row_timebase = timebase/float(ds0.number_of_rows)
    number_of_rows = ds0.number_of_rows
    frame_time = row_timebase*number_of_rows

    # --- TTree ---
    pname, ptitle = 'chanall', 'pulse tree'
    pt = ROOT.TTree(pname,ptitle)
    ptfill = pt.Fill
    ptbranch = pt.Branch
    bp_ch              = np.zeros(1,dtype=np.intc)
    bp_ev              = np.zeros(1,dtype=np.intc)
    bp_good            = np.zeros(1,dtype=bool)
    bp_grt             = np.zeros(1,dtype=np.intc)
    bp_prim            = np.zeros(1,dtype=bool)
    bp_shift1          = np.zeros(1,dtype=bool)
    bp_filt_phase      = np.zeros(1,dtype=np.float64)
    bp_filt_value      = np.zeros(1,dtype=np.float64)
    bp_filt_value_dc   = np.zeros(1,dtype=np.float64)
    bp_filt_value_phc  = np.zeros(1,dtype=np.float64)
    bp_filt_value_tdc  = np.zeros(1,dtype=np.float64)
    bp_energy          = np.zeros(1,dtype=np.float64)
    bp_energy_dc       = np.zeros(1,dtype=np.float64)
    bp_energy_phc      = np.zeros(1,dtype=np.float64)
    bp_min_value       = np.zeros(1,dtype=np.float64)
    bp_peak_index      = np.zeros(1,dtype=np.intc)
    bp_peak_time       = np.zeros(1,dtype=np.float64)
    bp_peak_value      = np.zeros(1,dtype=np.float64)
    bp_postpeak_deriv  = np.zeros(1,dtype=np.float64)
    bp_pretrig_mean    = np.zeros(1,dtype=np.float64)
    bp_pretrig_rms     = np.zeros(1,dtype=np.float64)
    bp_promptness      = np.zeros(1,dtype=np.float64)
    bp_pulse_average   = np.zeros(1,dtype=np.float64)
    bp_pulse_rms       = np.zeros(1,dtype=np.float64)
    bp_rise_time       = np.zeros(1,dtype=np.float64)
    bp_timestamp       = np.zeros(1,dtype=np.float64)
    bp_rowcount        = np.zeros(1,dtype=np.float64)
    bp_rowmodp         = np.zeros(1,dtype=np.float64)
    # --- cut paramaters ---
    bp_peak_region_max  = np.zeros(1, dtype=np.float64)                        
    bp_peak_region_mean = np.zeros(1, dtype=np.float64)
    bp_pre_region_sum   = np.zeros(1, dtype=np.float64)
    bp_jbr_region_sum   = np.zeros(1, dtype=np.float64)
    if EXTTRIG:
        bp_rows_after_last_external_trigger_nrp = np.zeros(1, dtype=np.int64)
        bp_rows_until_next_external_trigger_nrp = np.zeros(1, dtype=np.int64)
        bp_beamOn                      = np.zeros(1, dtype=bool)
    if GRTRIG:
        bp_sec_pr_maxmean = np.zeros(1, dtype=np.float64)            
        bp_sec_pr_meanmean = np.zeros(1, dtype=np.float64)            
        bp_sec_pr_mean = np.zeros(1, dtype=np.float64)            
        bp_sec_enemean = np.zeros(1, dtype=np.float64)            
        bp_sec_enemax = np.zeros(1, dtype=np.float64)            

    # --- tree branch address
    ptbranch('ch',             bp_ch,              'channel number/I')
    ptbranch('ev',             bp_ev,              'event/I')
    ptbranch('good',           bp_good,            'good/O')
    ptbranch('primary',        bp_prim,            'primary/O')
    ptbranch('grouptrigger',   bp_grt,             'grouptrigger/I')
    ptbranch('shift1',         bp_shift1,          'shift1/O')
    ptbranch('filt_phase',     bp_filt_phase,      'filt_phase/D')
    ptbranch('filt_value',     bp_filt_value,      'filt_value/D')
    ptbranch('filt_value_dc',  bp_filt_value_dc,   'filt_value_dc/D')
    ptbranch('filt_value_phc', bp_filt_value_phc,  'filt_value_phc/D')
    ptbranch('filt_value_tdc', bp_filt_value_tdc,  'filt_value_tdc/D')
    ptbranch('energy',         bp_energy,          'energy/D')
    ptbranch('energy_dc',      bp_energy_dc,       'energy_dc/D')
    ptbranch('energy_phc',     bp_energy_phc,      'energy_phc/D')
    ptbranch('min_value',      bp_min_value,       'min_value/D')
    ptbranch('peak_index',     bp_peak_index,      'peak_index/I')
    ptbranch('peak_time',      bp_peak_time,       'peak_time/D')
    ptbranch('peak_value',     bp_peak_value,      'peak_value/D')
    ptbranch('postpeak_deriv', bp_postpeak_deriv,  'postpeak_deriv/D')
    ptbranch('pretrig_mean',   bp_pretrig_mean,    'pretrig_mean/D')
    ptbranch('pretrig_rms',    bp_pretrig_rms,     'pretrig_rms/D')
    ptbranch('promptness',     bp_promptness,      'promptness/D')
    ptbranch('pulse_average',  bp_pulse_average,   'pulse_average/D')
    ptbranch('pulse_rms',      bp_pulse_rms,       'pulse_rms/D')
    ptbranch('rise_time',      bp_rise_time,       'rise_time/D')
    ptbranch('timestamp',      bp_timestamp,       'timestamp/D')
    ptbranch('rowcount',       bp_rowcount,        'rowcount/D')
    ptbranch('rowmodp',        bp_rowmodp,         'rowmodifiedp/D')
    # --- cut paramaters 
    ptbranch('peak_region_max',  bp_peak_region_max,  'peak_region_max/D')
    ptbranch('peak_region_mean', bp_peak_region_mean, 'peak_region_mean/D')
    ptbranch('pre_region_sum',   bp_pre_region_sum,   'pre_region_sum/D')
    ptbranch('jbr_region_sum',   bp_jbr_region_sum,   'jbr_region_sum/D')
    if EXTTRIG:
        ptbranch('row_after_extrig_nrp',bp_rows_after_last_external_trigger_nrp, 'row_after_extrig_nrp/L')
        ptbranch('row_next_extrig_nrp', bp_rows_until_next_external_trigger_nrp, 'row_next_extrig_nrp/L')
        ptbranch('beam',            bp_beamOn,                      'beam/O')
    if GRTRIG:
        ptbranch('sec_pr_maxmean',  bp_sec_pr_maxmean,  'sec_pr_maxmean/D')
        ptbranch('sec_pr_meanmean', bp_sec_pr_meanmean, 'sec_pr_meanmean/D')
        ptbranch('sec_pr_mean',     bp_sec_pr_mean,     'sec_pr_mean/D')
        ptbranch('sec_enemean',     bp_sec_enemean,     'sec_enemean/D')
        ptbranch('sec_enemax',      bp_sec_enemax,      'sec_enemax/D')
        
    # loop start
    for ds in data:
        dschan=ds.channum
        print "channel %d start.... for %.3f (sec)"%(dschan, (time.time() - start))
        # ---------------------------------------------
        bc_ch[0] = dschan
        bc_col[0] = ds.column_number
        bc_row[0] = ds.row_number
        bc_ncols[0] = ds.number_of_columns
        bc_nrows[0] = ds.number_of_rows
        bc_npulses[0] = ds.nPulses
        bc_nsamples[0] = ds.nSamples
        bc_npresamples[0] = ds.nPresamples
        bc_row_timebase[0] = ds.timebase/float(ds.number_of_rows)
        bc_timebase[0] = ds.timebase
        bc_timestamp_offset[0] = ds.timestamp_offset
        ctfill()
        # ---------------------------------------------            
        # --- numpy array objects for fast calculation
        np_good          =np.array(ds.good())
        np_prim          =np.zeros(shape=np_good.shape,dtype=bool)
        np_grt           =np.zeros(shape=np_good.shape,dtype=np.intc)
        np_shift1        =np.array(ds.p_shift1)
        np_filt_phase    =np.array(ds.p_filt_phase)
        np_filt_value    =np.array(ds.p_filt_value)
        np_filt_value_dc =np.array(ds.p_filt_value_dc)
        np_filt_value_phc=np.array(ds.p_filt_value_phc)
        np_filt_value_tdc=np.array(ds.p_filt_value_tdc)
        np_energy        =np.zeros(shape=np_filt_value.shape,dtype=np.float64)
        np_energy_dc     =np.zeros(shape=np_filt_value_dc.shape,dtype=np.float64)
        np_energy_phc    =np.zeros(shape=np_filt_value_phc.shape,dtype=np.float64)
        np_energy        = np.array(ds.p_energy)
        attr="p_filt_value_dc"
        if ds.calibration.has_key(attr):
            np_energy_dc = ds.calibration[attr].ph2energy(np_filt_value_dc)
        attr="p_filt_value_phc"
        if ds.calibration.has_key(attr):
            np_energy_phc = ds.calibration[attr].ph2energy(np_filt_value_phc)
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
        row_adjust = -1.*(np_filt_phase+np_shift1)
        np_rowmodp       =np_rowcount+(row_adjust*util.NUM_ROWS)-util.GLOBAL_PT_OFFSET
        # --- cut paramaters 
        np_peak_region_max =np.array(ds.p_peak_region_max)
        np_peak_region_mean=np.array(ds.p_peak_region_mean)
        np_pre_region_sum  =np.array(ds.p_pre_region_sum)
        np_jbr_region_sum  =np.array(ds.p_jbr_region_sum)
        # --- group trigger ----
        if GRTRIG:
            np_grt            =np.array(ds.p_grouptrig)
            prim_ind = np.where(np_grt==-1)[0]
            np_prim[prim_ind]=True
            np_sec_pr_maxmean =np.array(ds.sec_pr_maxmean)
            np_sec_pr_meanmean=np.array(ds.sec_pr_meanmean)
            np_sec_pr_mean    =np.array(ds.sec_pr_mean)
            np_sec_enemean    =np.array(ds.sec_enemean)
            np_sec_enemax     =np.array(ds.sec_enemax)
        # --- external trigger ----
        if EXTTRIG:
            np_rows_after_last_external_trigger_nrp=np.array(ds.rows_after_last_external_trigger_nrp)
            np_rows_until_next_external_trigger_nrp=np.array(ds.rows_until_next_external_trigger_nrp)
            np_beamOn                              =np.array(ds.beamflag)
        
        for i in xrange(ds.nPulses):
            bp_ch[0]              = dschan
            bp_ev[0]              = i
            bp_good[0]            = np_good[i]          
            bp_grt[0]             = np_grt[i]                     
            bp_prim[0]            = np_prim[i]                    
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
            bp_rowmodp[0]         = np_rowmodp[i]                 
            bp_peak_region_max[0] = np_peak_region_max[i]
            bp_peak_region_mean[0]= np_peak_region_mean[i]
            bp_pre_region_sum[0]  = np_pre_region_sum[i]
            bp_jbr_region_sum[0]  = np_jbr_region_sum[i]
            if GRTRIG:
                bp_sec_pr_maxmean[0]  = np_sec_pr_maxmean[i]
                bp_sec_pr_meanmean[0] = np_sec_pr_meanmean[i]
                bp_sec_pr_mean[0]     = np_sec_pr_mean[i]
                bp_sec_enemean[0]     = np_sec_enemean[i]
                bp_sec_enemax[0]      = np_sec_enemax[i]
            if EXTTRIG:
                bp_rows_after_last_external_trigger_nrp[0] = np_rows_after_last_external_trigger_nrp[i]
                bp_rows_until_next_external_trigger_nrp[0] = np_rows_until_next_external_trigger_nrp[i]
                bp_beamOn[0]                      = np_beamOn[i]
            ptfill()
            
    pt.Write()# event loop end
    ct.Write()# ds loop end

    f.Close()
    print "End of dump to ROOT for %.3f sec"%((time.time() - start))
    print "-------------------------------------------------"
    
