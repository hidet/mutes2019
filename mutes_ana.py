""" mutes_ana.py describes MUTES class

History: 
2018-07-04 ; ver 1.0; created from KHE.py, E62 run 
2018-07-05 ; ver 1.1; minor update by S.Y. 
2018-07-06 ; ver 1.2; TTree format changed from ch branch to single branch, S.Y.
2018-08-17 ; ver 1.3; file divided to mutes_ext,_group,_dump, T.H.
2019-01-31 ; ver 1.4; drastically simplified, especially removed plot functions by a bad boy HT
2019-05-09 ; ver 1.5; HT minor change
2019-05-23 ; ver 1.6; HT bug fixed
2019-06-04 ; ver 1.7; HT added row decimals

"""

__version__ = '1.7'

import mass
import math
import numpy as np
import pylab as plt

params = {'xtick.labelsize': 10, # x ticks
          'ytick.labelsize': 10, # y ticks
          'legend.fontsize': 8
                    }
plt.rcParams['font.family'] = 'serif' # 
plt.rcParams.update(params)

import os, sys
import h5py
import time
from collections import OrderedDict
import mutes_util as util
import mutes_ext as ext
import mutes_group as grp
import mutes_dump as dump
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
import matplotlib.cm as cm
import pprint
import tesmap_forgrptrig as tesmap
#plt.ion()
util = reload(util)
ext = reload(ext)
grp = reload(grp)
dump = reload(dump)
tesmap = reload(tesmap)


class MUTES():
    def __init__(self,pulse_runnums, noise_runnum, maxchans,
                 calibration_runnum, calibration_noisenum, badchans,
                 DATADIR, DELETE, GRTINFO, COLUMN_INFO,
                 catecut=None, target="Mn", cut_pre=0, cut_post=0,
                 use_new_filters=True):
        self.noise_runnum=noise_runnum
        self.calibration_runnum=calibration_runnum
        self.calibration_noisenum=calibration_noisenum
        self.maxchans=maxchans
        self.DATADIR = DATADIR
        self.GRTINFO = GRTINFO
        self.COLUMN_INFO = COLUMN_INFO
        self.badchans = badchans
        if isinstance(pulse_runnums,list)==False:
            pulse_runnums = (pulse_runnums,)
            self.pulse_runnums=tuple(pulse_runnums)
        self.pulse_runnums=pulse_runnums
        self.pulse_files,self.noise_files = util.get_multiple_file_lists(self.pulse_runnums,self.noise_runnum,
                                                                         self.maxchans,self.badchans,self.DATADIR)
        self.first_pulse_file = self.pulse_files[0][0]
        self.first_pulse_runnum=self.pulse_runnums[0]
        self.multiruns=""
        if len(self.pulse_runnums)>1:
            for run in self.pulse_runnums:
                self.multiruns += "_%s"%(run)
        self.usechans = util.get_usechans(self.first_pulse_runnum,self.noise_runnum,self.maxchans,self.badchans,self.DATADIR)

        self.cut_pre = cut_pre
        self.cut_post = cut_post
        self.use_new_filters=use_new_filters
        self.target = target
        self.catecut = catecut
        self.spill=None
        self.extall=True
        self.linefit_dict = {}

        self.trans = False
        self.calibration_hdf5_filename = None
        if not self.calibration_runnum is None:
            self.trans = True
            cal_pulse_files,cal_noise_files = util.get_file_lists(self.calibration_runnum,self.noise_runnum,
                                                                  self.maxchans,self.badchans,self.DATADIR)
            self.calibration_hdf5_filename = util.generate_hdf5_filename(cal_pulse_files[0],"_noi%04d"%self.calibration_noisenum+"_mass_2019") # assuming add="" in the name of a calibration run
       
        add=self.multiruns
        if self.use_new_filters==False:           add += "_old_filter"
        if self.trans:                            add += "_trans%d"%self.calibration_runnum
        if self.cut_pre>0 or self.cut_post>0:     add += "_pre%03d_post%03d"%(self.cut_pre,self.cut_post)
        if len(self.catecut)>0:
            if self.catecut.has_key('beam'):
                if self.catecut['beam']=='off':   add += "_spilloff"
                elif self.catecut['beam']=='on':  add += "_spillon"
                elif self.catecut['beam']=='None': self.catecut.pop('beam')
            if self.catecut.has_key('sprmc'):
                if self.catecut['sprmc']=='off':  add += "_sprmcoff"
                elif self.catecut['sprmc']=='on': add += "_sprmcon"
                elif self.catecut['sprmc']=='None': self.catecut.pop('sprmc')
            if self.catecut.has_key('jbrsc'):
                if self.catecut['jbrsc']=='off':  add += "_jbrscoff"
                elif self.catecut['jbrsc']=='on': add += "_jbrscon"
                elif self.catecut['jbrsc']=='None': self.catecut.pop('jbrsc')
            if self.catecut.has_key('prime'):
                if self.catecut['prime']=='off':  add += "_sec"
                elif self.catecut['prime']=='on': add += "_prime"
                elif self.catecut['prime']=='None': self.catecut.pop('prime')
            self.catecutname=str(catecut).replace(":","").replace(" ","").replace("'","").replace("{","").replace("}","").replace(",","_")
        if len(self.catecut)==0:
            self.catecut=None
            self.catecutname="nocatecut" 
                
        self.hdf5_filename       = util.generate_hdf5_filename(self.first_pulse_file,"_noi%04d"%self.noise_runnum+"_mass_2019"+add)
        self.hdf5_noisefilename  = util.generate_hdf5_filename(self.first_pulse_file,"_noi%04d"%self.noise_runnum+"_noise"+add)
        if DELETE: self.delete_hdf5_outputs()

        self.rootdir=dump.ROOTDIR+"/run%04d/"%self.first_pulse_runnum
        if os.path.isdir(self.rootdir)==False: os.makedirs(self.rootdir)
        self.root_filename = util.generate_user_root_filename(self.rootdir,"run%04d_noi%04d"%(self.first_pulse_runnum,self.noise_runnum)+"_mass_2019"+add)
            
        self.data = mass.TESGroup(self.pulse_files, self.noise_files, hdf5_filename=self.hdf5_filename, hdf5_noisefilename=self.hdf5_noisefilename)
        
    def delete_hdf5_outputs(self):
        try:
            os.remove(self.hdf5_filename)
        except OSError:
            pass
        try:
            os.remove(self.hdf5_noisefilename)
        except OSError:
            pass

    def define_spill(self,threshold=0.4):
        ds = self.data.first_good_dataset
        util.init_row_timebase(ds)
        print "define spills..."
        self.spill = ext.Spill(ds,threshold=threshold)

    def adjust_rowcount(self):
        for ds in self.data:
            # be careful, to adjust phase, you do filter the data first
            p_row_adj = -1.*(np.asarray(ds.p_shift1,dtype=np.float64)+np.asarray(ds.p_filt_phase,dtype=np.float64))*util.NUM_ROWS
            tmp=np.array([math.modf(p) for p in p_row_adj])
            setattr(ds,"p_rowp",np.asarray(ds.p_rowcount+tmp[:,1], dtype=np.int64)-util.GLOBAL_PT_OFFSET)# plus
            setattr(ds,"p_rown",np.asarray(ds.p_rowcount-tmp[:,1], dtype=np.int64)-util.GLOBAL_PT_OFFSET)# minus
            setattr(ds,"p_rowd",np.asarray(tmp[:,0], dtype=np.float32))# decimal
        
    def get_basic_cuts(self):
        pave_high=10000.
        peak_value_max = 40000.
        cuts = mass.core.controller.AnalysisControl(
            pulse_average=(10, pave_high),
            pretrigger_rms=(1, 50),
            peak_value=(1000, peak_value_max),
            #        postpeak_deriv=(0, 50),
            rise_time_ms=(0.07, 0.25),
            peak_time_ms=(0.1,  0.8),
        )
        return cuts
        

    def ana(self,forceNew=False,summaryNew=False,calibNew=False,
            exttrigNew=False,grptrigNew=False,bcutflag=True):

        if forceNew==True: calibNew=True
        # --- set good for bad channels ---
        #self.data.set_chan_good(self.data.why_chan_bad.keys())
        for ds in self.data.datasets:
            if ds.channum not in self.data.good_channels:
                self.data.set_chan_good(ds.channum)
                    
        # ------ DO NOT USE cython with record length cuts ------ HT 20190116
        use_cython = False
        print 'summarize_data', summaryNew
        self.data.summarize_data(cut_pre=self.cut_pre,cut_post=self.cut_post,
                                 use_cython=use_cython,forceNew=summaryNew,
                                 use_new_filters=self.use_new_filters)

        print 'basic cuts', bcutflag
        if bcutflag==True:
            bcut = self.get_basic_cuts()
            self.data.apply_cuts(bcut, forceNew=True)
        for ds in self.data:
            h5_good = ds.hdf5_group.require_dataset("good",(ds.nPulses,),dtype=np.int32)
            h5_good[:] = np.array(ds.good()).astype(int)
            setattr(ds,"p_goodflag",ds.hdf5_group["good"][()])
        
        # please fill categorical cuts first (e.g., beam "on:off" or sprmc "on:off")
        self.prime_analysis()
        self.jbrs_analysis(jbrsc_th=350,jbrsc_thn=-400)
        # external trigger data checking
        for pr in self.pulse_runnums:
            if (ext.check_external_trigger_data(pr,self.DATADIR))==False:
                self.extall=False
                print "Error: external trigger file is missing on run %d"%pr
        if self.extall:
            #self.define_spill(threshold=0.4)
            self.beamflag_analysis(forceNew=exttrigNew)
        self.group_trigger_peak_region_analysis(forceNew=grptrigNew,sprmc_th=10,sprmc_thn=-10)

        # basic analysis
        if self.trans==False:
            print 'basic analysis', forceNew
            self.data.avg_pulses_auto_masks(forceNew=forceNew,category=self.catecut)
            self.data.compute_filters(cut_pre=self.cut_pre,cut_post=self.cut_post,forceNew=forceNew,category=self.catecut)
            self.data.filter_data(forceNew=forceNew)
            self.data.drift_correct(forceNew=forceNew,category=self.catecut)
            self.data.phase_correct(forceNew=forceNew,category=self.catecut)
            self.mass_calibration_analysis(forceNew=calibNew,category=self.catecut)
            self.adjust_rowcount() # ds.p_rowp, ds.p_rown
            if self.extall: self.timing_analysis(forceNew=exttrigNew) # external trigger timing
            self.group_trigger_peak_region_analysis_filtered(forceNew=grptrigNew) # secondary peak region ene
        else:
            self.mass_analysis_transfer_calibration(forceNew,exttrigNew,grptrigNew)

            
    def mass_calibration_analysis(self, forceNew=False, category=None):
        self.calib_list = []
        if self.calibration_runnum is None:
            nextra=3
            if self.target == "CrCoCu":
                calib_list = ["CrKAlpha","CrKBeta","CoKAlpha","CoKBeta","CuKAlpha"]
            elif self.target == "CrCo":
                calib_list = ["CrKAlpha","CrKBeta","CoKAlpha","CoKBeta"]
            elif self.target == "Mn":             
                calib_list = ["MnKAlpha","MnKBeta"]
            elif self.target == "Fe":             
                calib_list = ["FeKAlpha","FeKBeta"]
            elif self.target == "Co57":# added for Co57 analysis, pls check data_TMU_XXXX.csv file
                calib_list = ["FeKAlpha","FeKBeta","Co57_14keV"]
                mass.STANDARD_FEATURES["Co57_14keV"]=14412.95# 14.4 keV
                nextra=0
            else:
                print "[Error] calib_list is out of range ", self.target
                return False
            self.calib_list = calib_list
            attrs=["p_filt_value_dc","p_filt_value_phc"]
            # NOTE: self.p_energy is overwritten with the last used attr
            for attr in attrs:
                print "......      in mass_calibration_analysis : attr = ", attr
                print self.target, self.calib_list
                self.data.calibrate(attr,calib_list,size_related_to_energy_resolution=100,nextra=nextra,
                                    forceNew=forceNew,category=category,vary_bg=True,vary_tail=True)
                
    def mass_analysis_transfer_calibration(self, forceNew=False, exttrigNew=False, grptrigNew=False):
        calh5name = self.calibration_hdf5_filename
        if not os.path.isfile(calh5name):
	    print "%s is not found... end the calibration transfer"%calh5name
            print "Please create %s OR set calrun as None in the csv file"%calh5name
            return False
        caldata = mass.TESGroupHDF5(calh5name)
        
        for ds in self.data:
            if (caldata.channel.has_key(ds.channum) and
                not caldata.why_chan_bad.has_key(ds.channum)):
                dscal = caldata.channel[ds.channum]
                if not ("calibration" in dscal.hdf5_group and "filters" in dscal.hdf5_group):
                    self.data.set_chan_bad(ds.channum, "caldata was bad")
                    continue
            else:
                self.data.set_chan_bad(ds.channum, "caldata was bad")
                continue
            ds.hdf5_group.require_group("filters")
            ds.hdf5_group["filters"].attrs["shorten"]=dscal.hdf5_group["filters"].attrs["shorten"]
            ds.hdf5_group["filters"].attrs["newfilter"]=dscal.hdf5_group["filters"].attrs["newfilter"]
            try:
                ds.hdf5_group["filters"]["filt_aterms"]=dscal.hdf5_group["filters"]["filt_aterms"][()]
            except:
                ds.hdf5_group["filters"]["filt_aterms"][:]=dscal.hdf5_group["filters"]["filt_aterms"][()]
            try:
                ds.hdf5_group["filters"]["filt_noconst"]=dscal.hdf5_group["filters"]["filt_noconst"][()]
            except:
                ds.hdf5_group["filters"]["filt_noconst"][:]=dscal.hdf5_group["filters"]["filt_noconst"][()]
            ds.hdf5_group["average_pulse"][:]=dscal.hdf5_group["average_pulse"][()]
            ds._MicrocalDataSet__setup_vectors()
            ds.filter_data(forceNew=forceNew)
           
            ds.p_filt_value_dc.attrs["type"]                = "ptmean_gain"
            ds.p_filt_value_dc.attrs["median_pretrig_mean"] = dscal.p_filt_value_dc.attrs["median_pretrig_mean"]
            ds.p_filt_value_dc.attrs["slope"]               = dscal.p_filt_value_dc.attrs["slope"]
            ds._apply_drift_correction()

            hdf5_cal_group = dscal.hdf5_group["calibration"]
            for k in hdf5_cal_group.keys():
                ds.calibration[k] = mass.EnergyCalibration.load_from_hdf5(hdf5_cal_group, k)

        self.data.phase_correct(forceNew=forceNew,category=self.catecut)
        #self.data.convert_to_energy("p_filt_value_dc", calname="p_filt_value_dc")
        #self.data.convert_to_energy("p_filt_value_dc", calname="p_filt_value_phc")
        self.data.convert_to_energy("p_filt_value_phc", calname="p_filt_value_phc")
        self.adjust_rowcount() # ds.p_rowp, ds.p_rown
        if self.extall: self.timing_analysis(forceNew=exttrigNew) # external trigger timing
        self.group_trigger_peak_region_analysis_filtered(forceNew=grptrigNew) # secondary peak region ene


    def beamflag_analysis(self,forceNew=False):
        self.data.register_categorical_cut_field("beam",["on","off"])
        for ds in self.data:
            util.init_row_timebase(ds)
            if hasattr(ds,"external_trigger_rowcount"):
                #if self.spill is not None: ext.define_beam_timing(ds,self.spill,forceNew=forceNew)
                ext.define_beam_timing(ds,np.asarray(ds.external_trigger_rowcount[:],dtype=np.int64),forceNew=forceNew)
        
    def timing_analysis(self,forceNew=False):
        self.data.register_categorical_cut_field("p_dtflag",["on","off"])
        print "timing analysis starts...", forceNew
        for ds in self.data:
            util.init_row_timebase(ds)
            if hasattr(ds,"external_trigger_rowcount"):
                ext.calc_external_trigger_timing(ds,np.asarray(ds.external_trigger_rowcount[:],dtype=np.int64),forceNew=forceNew)
                #if self.spill is not None: ext.calc_spill_timing(ds,np.array(self.spill.spill_start_rowcount,dtype=np.int64),forceNew=forceNew)

    def prime_analysis(self):
        print "primary analysis starts... (always update)"
        self.data.register_categorical_cut_field("prime",["on","off"])
        for ds in self.data:
            grti = np.array(ds.p_grouptrig)# group trig ch info
            prima = np.ones(shape=grti.shape,dtype=bool)
            sec_ind = np.where(grti!=-1)[0]# -1 is primary event, other number is a secondary one
            prima[sec_ind]=False
            ds.cuts.cut_categorical("prime", {"on": prima,"off": ~prima})
            h5_prime = ds.hdf5_group.require_dataset("prime",(ds.nPulses,),dtype=np.int32)
            h5_prime[:] = prima.astype(int)
            setattr(ds,"p_prime",ds.hdf5_group["prime"][()])
            
    def jbrs_analysis(self,jbrsc_th=350,jbrsc_thn=-400):
        print "jbr analysis starts... (always update)"
        self.data.register_categorical_cut_field("jbrsc",["on","off"])
        for ds in self.data:
            jbrscp = np.array(ds.p_jbr_region_sum)<jbrsc_th
            jbrscn = np.array(ds.p_jbr_region_sum)>jbrsc_thn
            jbrsc = np.logical_and(jbrscp,jbrscn)
            ds.cuts.cut_categorical("jbrsc", {"on": jbrsc,"off": ~jbrsc})
            h5_jbrsc = ds.hdf5_group.require_dataset("jbrsc",(ds.nPulses,),dtype=np.int32)
            h5_jbrsc[:] = jbrsc.astype(int)
            setattr(ds,"p_jbrsc",ds.hdf5_group["jbrsc"][()])

    def group_trigger_peak_region_analysis(self,forceNew=False,sprmc_th=10,sprmc_thn=-10):
        """
        sprmc_th : sec_pr_mean < 37 is a preliminary value derived by run 269 (X-ray and Beam), run 136,137,138, and 139 (beam only)
        10 could be a good starting point after removed MUX neighbor cross talk, and added the limit of negative side
        """
        print "group trigger analysis starts...", forceNew
        self.data.register_categorical_cut_field("sprmc",["on","off"])
        flag=grp.calc_group_trigger_params(self.data,self.GRTINFO,forceNew=forceNew)
        if flag==False:
            "no group trigger data were set"
            return
        for ds in self.data:
            sprmcp = np.array(ds.sec_pr_mean)<sprmc_th# secondary peak region mean cut
            sprmcn = np.array(ds.sec_pr_mean)>sprmc_thn# secondary peak region mean cut negative side
            sprmc = np.logical_and(sprmcp,sprmcn)
            ds.cuts.cut_categorical("sprmc", {"on": sprmc,"off": ~sprmc})
            h5_sprmc = ds.hdf5_group.require_dataset("sprmc",(ds.nPulses,),dtype=np.int32)
            h5_sprmc[:] = sprmc.astype(int)# becareful: the categorical cut means logical_and of good and sprmc
            setattr(ds,"p_sprmc",ds.hdf5_group["sprmc"][()])
            
    def group_trigger_peak_region_analysis_filtered(self,forceNew=False):
        print "group trigger analysis filtered starts...", forceNew
        grp.calc_group_trigger_params_filtered(self.data,self.GRTINFO,forceNew=forceNew)

    def dump_ROOT(self, EXTTRIG=False, GRTRIG=False, dumppulse=False):
        dump.dump_ROOT(self.data,self.root_filename,EXTTRIG,GRTRIG,dumppulse)

