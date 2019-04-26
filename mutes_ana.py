""" mutes_ana.py describes KHE class used for E62 analysis 

History: 
2018-07-04 ; ver 1.0; created from KHE.py, E62 run 
2018-07-05 ; ver 1.1; minor update by S.Y. 
2018-07-06 ; ver 1.2; TTree format changed from ch branch to single branch, S.Y.
2018-08-17 ; ver 1.3; file divided to mutes_ext,_group,_dump, T.H.
2019-01-31 ; ver 1.4; drastically simplified, especially removed plot functions by a bad boy HT

"""

__version__ = '1.4'

import mass
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


class KHE():
    def __init__(self,pulse_runnums, noise_runnum, maxchans,
                 calibration_runnum, calibration_noisenum, badchans,
                 DATADIR, DELETE, GRTINFO, COLUMN_INFO,
                 hdf5optname=None, catecut=None, target="Mn", cut_pre=0, cut_post=0):
        self.pulse_runnums=pulse_runnums
        self.noise_runnum=noise_runnum
        self.calibration_runnum=calibration_runnum
        self.calibration_noisenum=calibration_noisenum
        self.linefit_dict = {}
        self.DATADIR = DATADIR
        self.GRTINFO = GRTINFO
        self.COLUMN_INFO = COLUMN_INFO
        self.cut_pre = cut_pre
        self.cut_post = cut_post
        self.badchans = badchans
        self.multiruns=""
        if isinstance(pulse_runnums,list):
            pulse_files,noise_files = util.get_multiple_file_lists(self.pulse_runnums,self.noise_runnum,
                                                                   maxchans,self.badchans,self.DATADIR)
            first_pulse_file = pulse_files[0][0]
            self.first_pulse_runnum=self.pulse_runnums[0]
            self.MULTI=True
            for run in self.pulse_runnums:
                self.multiruns += "_%s"%(run)
        else:
            pulse_files,noise_files = util.get_file_lists(self.pulse_runnums,self.noise_runnum,
                                                          maxchans,self.badchans,self.DATADIR)
            first_pulse_file = pulse_files[0]
            self.first_pulse_runnum=self.pulse_runnums
            self.MULTI=False
        
        self.usechans = util.get_usechans(self.first_pulse_runnum,self.noise_runnum,maxchans,self.badchans,self.DATADIR)
        self.target = target
        self.catecut = catecut
        if self.catecut is None: 
            self.catecutname = "nocatecut" 
        else:
            self.catecutname=str(catecut).replace(":","").replace(" ","").replace("'","").replace("{","").replace("}","").replace(",","_")
        self.bonly = False
        if not self.calibration_runnum is None:
            self.bonly = True
            cal_pulse_files,cal_noise_files = util.get_file_lists(self.calibration_runnum,self.noise_runnum,
                                                                  maxchans,self.badchans,self.DATADIR)
            self.calibration_hdf5_filename = util.generate_hdf5_filename(cal_pulse_files[0],"_noi%04d"%self.calibration_noisenum+"_mass_2019")
            # this hdf5 file name should be matched, something catecuts included (_pre,_post,spillon/off,sprmcon/off,etc...)
        else:
            self.calibration_hdf5_filename = None

        self.rootdir=dump.ROOTDIR+"/run%04d/"%self.first_pulse_runnum
        if os.path.isdir(self.rootdir)==False: os.makedirs(self.rootdir)
        if hdf5optname is None:
            add=self.multiruns
            if self.bonly:
                add += "_trans%d"%self.calibration_runnum
            if self.cut_pre>0 or self.cut_post>0:
                add += str("_pre%03d_post%03d" %(self.cut_pre,self.cut_post))
            if self.catecut is not None and self.catecut.has_key('beam'):
                if self.catecut['beam']=='off': add = add + "_spilloff"
                elif self.catecut['beam']=='on': add = add + "_spillon"
            if self.catecut is not None and self.catecut.has_key('sprmc'):
                if self.catecut['sprmc']=='off': add = add + "_sprmcoff"
                elif self.catecut['sprmc']=='on': add = add + "_sprmcon"
            if self.catecut is not None and self.catecut.has_key('jbrsc'):
                if self.catecut['jbrsc']=='off': add = add + "_jbrscoff"
                elif self.catecut['jbrsc']=='on': add = add + "_jbrscon"
            if self.catecut is not None and self.catecut.has_key('prime'):
                if self.catecut['prime']=='off': add = add + "_sec"
                elif self.catecut['prime']=='on': add = add + "_prime"
            self.hdf5_filename       = util.generate_hdf5_filename(first_pulse_file,"_noi%04d"%self.noise_runnum+"_mass_2019"+add)
            self.hdf5_noisefilename  = util.generate_hdf5_filename(first_pulse_file,"_noi%04d"%self.noise_runnum+"_noise"+add)
            self.root_filename       = util.generate_user_root_filename(self.rootdir,"run%04d_noi%04d"%(self.first_pulse_runnum,self.noise_runnum)+"_mass_2019"+add)
        else:
            self.hdf5_filename       = util.generate_hdf5_filename(first_pulse_file,"_noi%04d"%self.noise_runnum+"_mass_2019"+str(hdf5optname))
            self.hdf5_noisefilename  = util.generate_hdf5_filename(first_pulse_file,"_noi%04d"%self.noise_runnum+"_noise"+str(hdf5optname))
            self.root_filename       = util.generate_user_root_filename(self.rootdir,"run%04d_noi%04d"%(self.first_pulse_runnum,self.noise_runnum)+"_mass_2019"+str(hdf5optname))

        if DELETE: self.delete_hdf5_outputs()
        self.data = mass.TESGroup(pulse_files, noise_files, hdf5_filename=self.hdf5_filename, hdf5_noisefilename=self.hdf5_noisefilename)
        
    def delete_hdf5_outputs(self):
        try:
            os.remove(self.hdf5_filename)
        except OSError:
            pass
        try:
            os.remove(self.hdf5_noisefilename)
        except OSError:
            pass

    def get_basic_cuts(self):
        #def get_basic_cuts(bonly=True):
        print "bonly = ", self.bonly
        pave_high=10000.
        peak_value_max = 40000.
        # cuts_for_beamonly = mass.core.controller.AnalysisControl(
        #     pulse_average=(0, pave_high),
        #  ) # this is just used for run 360, 412 for test, S.Y.
        cuts_for_beamonly = mass.core.controller.AnalysisControl(
            pulse_average=(10, pave_high),
            pretrigger_rms=(1, 50),
            peak_value=(1000, peak_value_max),
            #        postpeak_deriv=(0, 1000),
            rise_time_ms=(0.07, 0.25),
            peak_time_ms=(0.1,  0.8),     
        ) 
        cuts = mass.core.controller.AnalysisControl(
            pulse_average=(10, pave_high),
            pretrigger_rms=(1, 50),
            peak_value=(1000, peak_value_max),
            #        postpeak_deriv=(0, 50),
            rise_time_ms=(0.07, 0.25),
            peak_time_ms=(0.1,  0.8),
        )
        #    if self.bonly:
        #        return cuts_for_beamonly
        #    else:
        return cuts
        
    # simple analysis by H.Tatsuno
    def anahide(self,forceNew=False,summaryNew=False,calibNew=False,exttrigNew=False,grptrigNew=False,bcutflag=True):
        if forceNew==True:# if filter is changed, everything will change
            calibNew=True
        # --- set good for bad channels ---
        #self.data.set_chan_good(self.data.why_chan_bad.keys())
        for ds in self.data.datasets:
            if ds.channum not in self.data.good_channels:
                self.data.set_chan_good(ds.channum)
                    
        # need to always re-run this to populate some fields from hdf5
        #use_cython = True if self.cut_pre==0 and self.cut_post==0 else False
        use_cython = False
        # ------ DO NOT USE cython with record length cuts ------ HT 20190116
        print 'summarize_data', summaryNew
        self.data.summarize_data(cut_pre=self.cut_pre,cut_post=self.cut_post,use_cython=False,forceNew=summaryNew)

        print 'basic cuts', bcutflag
        if bcutflag==True:
            bcut = self.get_basic_cuts()
            self.data.apply_cuts(bcut, forceNew=True)
                
        for ds in self.data:
            h5_good = ds.hdf5_group.require_dataset("good",(ds.nPulses,),dtype=np.int32)
            h5_good[:] = np.array(ds.good()).astype(int)
            ds.goodflag = ds.hdf5_group["good"]
        
        # please fill categorical cuts first (e.g., beam "on:off" or sprmc "on:off")
        self.prime_analysis()
        self.jbrs_analysis(jbrsc_th=350,jbrsc_thn=-400)
        extall=True
        if self.MULTI:
            for pr in self.pulse_runnums:
                if (ext.check_external_trigger_data(pr))==False:
                    extall=False
                    print "Error: external trigger file is missing on run %d"%pr
        else:
            if (ext.check_external_trigger_data(self.pulse_runnums))==False:
                extall=False
                print "Error: external trigger file is missing on run %d"%self.pulse_runnums
        if extall:  self.timing_analysis(forceNew=exttrigNew)
        self.group_trigger_peak_region_analysis(forceNew=grptrigNew,sprmc_th=10,sprmc_thn=-10)
        # basic analysis
        if self.bonly==False:
            print 'basic analysis', forceNew
            self.data.avg_pulses_auto_masks(forceNew=forceNew,category=self.catecut)
            self.data.compute_filters(cut_pre=self.cut_pre,cut_post=self.cut_post,forceNew=forceNew,category=self.catecut)
            self.data.filter_data(forceNew=forceNew)
            self.data.drift_correct(forceNew=forceNew,category=self.catecut)
            self.data.phase_correct(forceNew=forceNew,category=self.catecut)
            self.mass_calibration_analysis(forceNew=calibNew,category=self.catecut)
        else:
            self.mass_analysis_transfer_calibration(forceNew=forceNew)

    def mass_calibration_analysis(self, forceNew=False, category=None):
        data = self.data
        self.calib_list = []
        if self.calibration_runnum is None:
            if self.target == "CrCoCu":
                calib_list = ["CrKAlpha","CrKBeta","CoKAlpha","CoKBeta","CuKAlpha"]
            elif self.target == "CrCo":
                calib_list = ["CrKAlpha","CrKBeta","CoKAlpha","CoKBeta"]
            elif self.target == "Mn":             
                calib_list = ["MnKAlpha","MnKBeta"]
            elif self.target == "Fe":             
                calib_list = ["FeKAlpha","FeKBeta"]
            else:
                print "[Error] calib_list is out of range ", self.target 
            self.calib_list = calib_list
            attrs=["p_filt_value_dc","p_filt_value_phc"]
            # NOTE: self.p_energy is overwritten with the last used attr
            for attr in attrs:
                print "......      in mass_calibration_analysis : attr = ", attr
                data.calibrate(attr,calib_list,size_related_to_energy_resolution=100,
                               forceNew=forceNew,category=category)

    def mass_analysis_transfer_calibration(self, forceNew=False):
        calh5name = self.calibration_hdf5_filename
        if not os.path.isfile(calh5name):
	    print "%s is not found..."%calh5name
            return
        self.caldata = mass.TESGroupHDF5(calh5name)
        for ds in self.data:
            if (self.caldata.channel.has_key(ds.channum) and
                not self.caldata.why_chan_bad.has_key(ds.channum)):
                dscal = self.caldata.channel[ds.channum]
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
                ds.hdf5_group["filters"]["filt_aterms"]=dscal.hdf5_group["filters"]["filt_aterms"].value
            except:
                ds.hdf5_group["filters"]["filt_aterms"][:]=dscal.hdf5_group["filters"]["filt_aterms"].value
            try:
                ds.hdf5_group["filters"]["filt_noconst"]=dscal.hdf5_group["filters"]["filt_noconst"].value
            except:
                ds.hdf5_group["filters"]["filt_noconst"][:]=dscal.hdf5_group["filters"]["filt_noconst"].value
            ds.hdf5_group["average_pulse"][:]=dscal.hdf5_group["average_pulse"].value
            ds._MicrocalDataSet__setup_vectors()
            ds.filter_data(forceNew=forceNew)
            
            ds.p_filt_value_dc.attrs["type"] = "ptmean_gain"
            ds.p_filt_value_dc.attrs["median_pretrig_mean"] = dscal.p_filt_value_dc.attrs["median_pretrig_mean"]
            ds.p_filt_value_dc.attrs["slope"] = dscal.p_filt_value_dc.attrs["slope"]
            ds._apply_drift_correction()
            
            hdf5_cal_group = dscal.hdf5_group["calibration"]
            for k in hdf5_cal_group.keys():
                ds.calibration[k] = mass.EnergyCalibration.load_from_hdf5(hdf5_cal_group, k)

        #self.data.convert_to_energy("p_filt_value_dc", calname="p_filt_value_dc")
        self.data.convert_to_energy("p_filt_value_dc", calname="p_filt_value_phc")
        

    def timing_analysis(self,forceNew=False):
        print "External trigger timing analysis starts...", forceNew
        self.data.register_categorical_cut_field("beam",["on","off"])
        ds = self.data.first_good_dataset
        util.init_row_timebase(ds)
        self.external_trigger_rowcount = np.asarray(ds.external_trigger_rowcount[:], dtype=np.int64)
        self.external_trigger_rowcount_as_seconds = np.array(ds.external_trigger_rowcount_as_seconds[:])
        #print "define spills..."
        #self.spill = ext.Spill(ds)
        for ds in self.data:
            util.init_row_timebase(ds)
            #ext.calc_external_trigger_timing(ds,self.spill,forceNew=forceNew)
            ext.calc_external_trigger_timing(ds,self.external_trigger_rowcount,forceNew=forceNew)

    def prime_analysis(self):
        print "KHE prime analysis starts... (always update)"
        self.data.register_categorical_cut_field("prime",["on","off"])
        for ds in self.data:
            grti = np.array(ds.p_grouptrig)# group trig ch info
            prima = np.zeros(shape=grti.shape,dtype=bool)
            prim_ind = np.where(grti==-1)[0]# -1 is primary event, other number is a secondary one
            prima[prim_ind]=True
            ds.cuts.cut_categorical("prime", {"on": prima,"off": ~prima})
            h5_prime = ds.hdf5_group.require_dataset("prime",(ds.nPulses,),dtype=np.int32)
            h5_prime[:] = prima.astype(int)
            
    def jbrs_analysis(self,jbrsc_th=350,jbrsc_thn=-400):
        print "KHE jbr analysis starts... (always update)"
        self.data.register_categorical_cut_field("jbrsc",["on","off"])
        for ds in self.data:
            jbrscp = np.array(ds.p_jbr_region_sum)<jbrsc_th
            jbrscn = np.array(ds.p_jbr_region_sum)>jbrsc_thn
            jbrsc = np.logical_and(jbrscp,jbrscn)
            ds.cuts.cut_categorical("jbrsc", {"on": jbrsc,"off": ~jbrsc})
            h5_jbrsc = ds.hdf5_group.require_dataset("jbrsc",(ds.nPulses,),dtype=np.int32)
            h5_jbrsc[:] = jbrsc.astype(int)
            
    def group_trigger_peak_region_analysis(self,forceNew=False,sprmc_th=10,sprmc_thn=-10):
        """
        sprmc_th : sec_pr_mean < 37 is a preliminary value derived by run 269 (X-ray and Beam), run 136,137,138, and 139 (beam only)

        History 
        2018/06/XX, H.Tatsuno, first draft 
        2018/07/05, S.Yamada, updated to use 37 as a default threshoud. 
        2019/01/20, H.Tatsuno, 10 is a good starting point after removed MUX neighbor, and added the limit of negative side
        """        
        print "group trigger analysis starts...", forceNew
        if self.catecutname == "nocatecut":
            pass
        else:
            self.catecutname = "%s_th%03d" % (self.catecutname, int(sprmc_th))

        print "forceNew=", forceNew,", catecutname = ", self.catecutname
        self.data.register_categorical_cut_field("sprmc",["on","off"])
        flag=grp.calc_group_trigger_params(self.data,self.GRTINFO,forceNew=forceNew)
        if not flag:
            print "Warning: group trigger analysis seems to be strange."
            return
        for ds in self.data:
            sprmcp = np.array(ds.sec_pr_mean)<sprmc_th# secondary peak region mean cut
            sprmcn = np.array(ds.sec_pr_mean)>sprmc_thn# secondary peak region mean cut negative side
            sprmc = np.logical_and(sprmcp,sprmcn)
            ds.cuts.cut_categorical("sprmc", {"on": sprmc,"off": ~sprmc})
            # fill the bool values of "sprmc" to the HDF5 file
            h5_sprmc = ds.hdf5_group.require_dataset("sprmc",(ds.nPulses,),dtype=np.int32)
            h5_sprmc[:] = sprmc.astype(int)# becareful: the categorical cut means logical_and of good and sprmc

    def dump_ROOT_2019(self, EXTTRIG=False, GRTRIG=False, dumppulse=False):
        dump.dump_ROOT_2019(self.data,self.root_filename,EXTTRIG,GRTRIG,dumppulse)

    def get_output_base(self):
        dirname = "{}/output/run{:04d}_n{:04d}".format(self.DATADIR, self.first_pulse_runnum,self.noise_runnum)
        if self.catecutname == "nocatecut":
            pass 
        else:
            dirname += "_"
            dirname += self.catecutname
        if self.cut_pre > 0 or self.cut_post>0:
            dirname += str("/pre%03d_post%03d" %(self.cut_pre,self.cut_post))

        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        return dirname+"/"

    def savefigs(self,closeFigs=True):
        for i in plt.get_fignums():
            plt.figure(i)
            plt.savefig(self.get_output_base()+"{}.png".format(i))
            if closeFigs: plt.close(i)

    def savetxt_linefit_each(self,linename):
        self.linename_list = self.linefit_dict.keys()
        self.linename = linename
        if not self.linename in self.linename_list:
            print "linename:%s has not plotted yet."%(self.linename)
            sys.exit()
        self.each_csv = self.get_output_base()+"run{:04d}_n{:04d}".format(self.first_pulse_runnum,self.noise_runnum)+ "_%s.csv"%self.linename
        fout = open(self.each_csv,'w')
        fout.write("chan,on_resol[eV],on_tail[%],off_resol[eV],off_tail[%]\n")       
        for chan in self.linefit_dict[linename].keys():
            fout.write("%d,%.2f,%.0f,%.2f,%.0f\n"%(chan,self.linefit_dict[self.linename][chan][0],self.linefit_dict[self.linename][chan][1]*100,self.linefit_dict[self.linename][chan][2],self.linefit_dict[self.linename][chan][3]*100))
        fout.close()
        print "result saved as %s"%(self.each_csv)
	

        
