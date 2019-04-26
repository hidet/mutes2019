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
            self.first_pulse_runnum=pulse_runnums[0]
            self.MULTI=True
            for run in self.pulse_runnums:
                self.multiruns += "_%s"%(run)
        else:
            pulse_files,noise_files = util.get_file_lists(self.pulse_runnums,self.noise_runnum,
                                                          maxchans,self.badchans,self.DATADIR)
            first_pulse_file = pulse_files[0]
            self.first_pulse_runnum=pulse_runnums
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
        for pr in self.pulse_runnums:
            if (ext.check_external_trigger_data(pr))==False:
                extall=False
                print "Error: external trigger file is missing on run %d"%pr

        # merge the external trigger files
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
        #external_trigger_rowcount=util.get_beam_clocks('clock2',self.first_pulse_runnum)
        #external_trigger_rowcount2=util.get_beam_clocks('clock3',self.first_pulse_runnum)
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
        print "KHE group trigger analysis starts...", forceNew
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
	

        
#    def setcut(self, cut_category, cutMin, cutMax, energy_cut,
#               energyMin=6900., energyMax=6975.,pulse_ana=False):
#        """
#        This is for setting cut by category (e.g. energy)
#        """
#        cutarray = []
#        ene_cutarray = []
# 
#        print "setcut with %s (pulse_ana=%s)"%(cut_category,str(pulse_ana))
#        for ds in self.data:
#            if pulse_ana == True:
#                if ds.channum in self.prim_chan_list:
##                if ds.channum in self.prim_chan_list[0]:
#                    if energy_cut:
#                        energy_array = np.array(ds.p_energy)
#                        ene_min = energy_array > energyMin
#                        ene_max = energy_array < energyMax
#                        tmp_ene_cutarray = ene_max * ene_min
#                    else:
#                        tmp_ene_cutarray = True # it means all true
#                        # energy_array = np.array(ds.p_energy)
#                        # tmp_ene_cutarray = np.ones(shape=energy_array.shape,dtype=bool)
#                    if cut_category == "None":
#                        tmp_cut = True
#                        #tmp_cut = np.ones(shape=energy_array.shape,dtype=bool)
#                    elif cut_category == "sec_pr_mean":
#                        if  cutMin != None and cutMax != None:
#                            sec_pr_mean_array = np.array(ds.sec_pr_mean)
#                            cut_max = sec_pr_mean_array < cutMax
#                            cut_min = sec_pr_mean_array > cutMin
#                            tmp_cut = cut_max * cut_min
#                        elif cutMin == None:
#                            sec_pr_mean_array = np.array(ds.sec_pr_mean)
#                            cut_max = sec_pr_mean_array < cutMax
#                            tmp_cut = cut_max
#                        elif cutMax == None:
#                            sec_pr_mean_array = np.array(ds.sec_pr_mean)
#                            cut_min = sec_pr_mean_array > cutMin
#                            tmp_cut = cut_min
#                    
#                    tmp_cutarray = np.logical_and(tmp_cut,tmp_ene_cutarray)
#
#
#                else: continue
#            else:
#                
#                # if energy_cut:
#                #         energy_array = np.array(ds.p_energy)
#                #         ene_min = energy_array >= energyMin
#                #         ene_max = energy_array <= energyMax
#                #         tmp_ene_cutarray = ene_max * ene_min
#                # else:
#                #     ene_cutarray = True # it means all true
#            
#                if cut_category == "None":
#                    # tmp_cutarray = np.logical_and(tmp_ene_cutarray, True)
#                    tmp_cutarray = np.ones(shape=np.array(ds.p_energy).shape, dtype=bool)
#                elif cut_category == "sec_pr_mean":
#
#                    if cutMin != None and cutMax != None:
#                        sec_pr_mean_array = np.array(ds.sec_pr_mean)
#                        cut_max = sec_pr_mean_array < cutMax
#                        cut_min = sec_pr_mean_array > cutMin
#                        tmp_cutarray = cut_max * cut_min
#                    elif cutMin == None:
#                        sec_pr_mean_array = np.array(ds.sec_pr_mean)
#                        cut_max = sec_pr_mean_array < cutMax
#                        tmp_cutarray = cut_max
#                    elif cutMax == None:
#                        sec_pr_mean_array = np.array(ds.sec_pr_mean)
#                        cut_min = sec_pr_mean_array > cutMin
#                        tmp_cutarray = cut_min
#
#            cutarray.append(tmp_cutarray)
#            # ene_cutarray.append(tmp_ene_cutarray)
#        return cutarray, cut_category
#
#
#    
#            
#    def plot_all(self):# plot for w/o beam
#        print "printing figures..."
#        for linename in self.calib_list:
#            if 'KAlpha' in linename:
#                self.plot_linefit(linename)
#            if 'MnKBeta' in linename:
#                self.plot_linefit(linename)
#        self.plot_overall_hist()
#
#
#    def plot_all_beam(self, Docut=False, grptrigNew=False):# plot for with beam and spill
#        print "printing figures..."
#        self.spill.plot_all()
#        for linename in self.calib_list:
#            if 'KAlpha' in linename:
#                self.plot_linefit_beam(linename, grptrigNew = grptrigNew)
#        self.plot_check_beam_cut()
#        fname = self.get_output_base()+"2dhist_row_merged.hdf5"
##        ext.plot_and_save_hist2d_merged(data=self.data,fname=fname,tlo=-150, thi=150, tbin=1)
#        if Docut:
#            fname = self.get_output_base()+"2dhist_row_merged_with_cut.hdf5"
# #           ext.plot_and_save_hist2d_merged_with_cut(data=self.data,fname=fname, cutarray=self.cutarray,tlo=-150, thi=150, tbin=1)
#        self.plot_overall_hist_beam()
#
#    def plot_each(self):
#        for linename in self.calib_list:
#            if 'KAlpha' in linename:
#                self.plot_linefit_each(linename)
#
#    def plot_each_beam(self):
#        for linename in self.calib_list:
#            if 'KAlpha' in linename:
#                self.plot_linefit_each_beam(linename)
#
#    def plot_each_chan_resol_and_hist(self):
#        for linename in self.calib_list:
#            if 'KAlpha' in linename:
#                self.plot_each_resol_hist(linename)
#
#    def plot_each_chan_resol_and_hist_beam(self):
#        for linename in self.calib_list:
#            if 'KAlpha' in linename:
#                self.plot_each_resol_hist_beam(linename)
#
#    def savetxt_each(self):
#        for linename in self.calib_list:
#            if 'KAlpha' in linename:
#                self.savetxt_linefit_each(linename)
#
#    def plot_check_beam_cut(self):
#        if not self.spill.have_enough_spills():
#            return
#        ds = self.data.first_good_dataset
#        t = ds.external_trigger_rowcount
#        # lets plot the first n spills
#        n_spill = 40
#        n = self.spill.spill_end_inds[n_spill]
#        plt.figure()
#        median_energy = np.median(ds.p_energy[ds.good(beam="on")])
#        plt.plot(t[:n],median_energy*(t[:n]>0),".",label="external trigger")
#        plt.plot([t[i] for i in self.spill.spill_start_inds[:n_spill]], median_energy*np.ones(n_spill),"r.",label="start")
#        plt.plot([t[i] for i in self.spill.spill_end_inds[:n_spill]], median_energy*np.ones(n_spill),"m.",label="end")
#        plt.plot(ds.p_rowcount[ds.good(beam="on")],ds.p_energy[ds.good(beam="on")],".")
#        plt.plot(ds.p_rowcount[ds.good(beam="off")],ds.p_energy[ds.good(beam="off")],".")
#        plt.xlabel("time (s)")
#        plt.ylabel("spill status")
#        plt.legend(loc="best")
#        plt.title("{}, chan {:d}\nzoom in to convince yourself the beam on/off cut is right".format(self.spill.runstr,self.spill.ds.channum))
#
#    def plot_linefit(self,linename):
#        elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
#        edges = np.arange(elo,ehi,1)
#        counts, centers = util.hist_data(self.data, edges, category={})
#        fitter = util.linefit(counts, edges, linename)
#        plt.figure()
#        fitter.plot(axis=plt.gca(),label=False, color="r", ph_units="eV")
#        plt.xlabel("energy (eV)")
#        plt.ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#        res = fitter.last_fit_params_dict["resolution"][0]
#        tail_frac = fitter.last_fit_params_dict["tail_frac"][0]
#        plt.legend(["no beam","no beam: {:.2f} eV, {:.0f}% tail".format(res, tail_frac*100)],loc="best")
#        plt.title("{}, {} good chans\n{}".format(self.spill.runstr,self.data.n_good_channels(),linename))
#        
#
#    def plot_linefit_beam(self,linename, grptrigNew = False, detail=True):
#        elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
#        edges = np.arange(elo,ehi,1)
#        countson, centers = util.hist_data(self.data, edges, category={"beam":"on"})
#        countsoff, centers = util.hist_data(self.data, edges, category={"beam":"off"})
#
#        fitteron = util.linefit(countson, edges, linename)
#        fitteroff = util.linefit(countsoff, edges, linename)
#
#        if grptrigNew:
#            countsoncut, centers = util.hist_data(self.data, edges, category={"beam":"on", "sprmc":"on"})
#            fitteroncut = util.linefit(countsoncut, edges, linename)
#
#        plt.figure()
#        fitteron.plot(axis=plt.gca(),label=False, color="r", ph_units="eV")
#
#        if grptrigNew:
#            fitteroncut.plot(axis=plt.gca(),label=False, color="g", ph_units="eV")
#
#        fitteroff.plot(axis=plt.gca(),label=False, color="b", ph_units="eV")
#
#        plt.xlabel("energy (eV)")
#        plt.ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#
#        reson = fitteron.last_fit_params_dict["resolution"][0]
#        tail_frac_on = fitteron.last_fit_params_dict["tail_frac"][0]
#
#        resoff = fitteroff.last_fit_params_dict["resolution"][0]
#        tail_frac_off = fitteroff.last_fit_params_dict["tail_frac"][0]
#
#        plt.legend(["beam on","on: {:.2f} eV, {:.0f}% tail".format(reson, tail_frac_on*100),
#                    "beam off","off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)],loc="best")
#
#        if grptrigNew:
#            resoncut = fitteroncut.last_fit_params_dict["resolution"][0]
#            tail_frac_oncut = fitteroncut.last_fit_params_dict["tail_frac"][0]
#
#            plt.legend(["beam on","on: {:.2f} eV, {:.0f}% tail".format(reson, tail_frac_on*100),
#                        "beam on cut, sprmc_th < " + str(self.sprmc_th),"on: {:.2f} eV, {:.0f}% tail".format(resoncut, tail_frac_oncut*100),
#                        "beam off","off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)],loc="best")
#            fout = open(self.get_output_base()+linename+"_resol.csv",'w')
#            fout.write("%d,%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n"%(self.first_pulse_runnum,self.cut_pre,self.cut_post,self.data.n_good_channels(),reson,tail_frac_on,resoncut,tail_frac_oncut,resoff,tail_frac_off))
#            fout.close()
#
#        plt.title("{}, {} good chans\n{}".format(self.spill.runstr,self.data.n_good_channels(),linename))
#
#        if detail:
#            if grptrigNew:
#                plt.figure()
#                fitteroncut.plot(axis=plt.gca(),label=False, color="g", ph_units="eV")
#                fitteroff.plot(axis=plt.gca(),label=False, color="b", ph_units="eV")
#
#                plt.xlabel("energy (eV)")
#                plt.ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#                #plt.legend(["beam on cut, sprmc_th < " + str(self.sprmc_th),"on: {:.2f} eV, {:.0f}% tail".format(resoncut, tail_frac_oncut*100),
##                        "beam off","off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)],loc="best")
#                plt.title("{}, {} good chans\n{}".format(self.spill.runstr,self.data.n_good_channels(),linename))
#
#                plt.figure()
#                fitteroncut.plot(axis=plt.gca(),label=False, color="g", ph_units="eV")
#
#                plt.xlabel("energy (eV)")
#                plt.ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#                #plt.legend(["beam on cut, sprmc_th < " + str(self.sprmc_th),"on: {:.2f} eV, {:.0f}% tail".format(resoncut, tail_frac_oncut*100),
##                        "beam off","off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)],loc="best")
#                plt.title("{}, {} good chans\n{}".format(self.spill.runstr,self.data.n_good_channels(),linename))
#
#
#                plt.figure()
#                fitteroff.plot(axis=plt.gca(),label=False, color="b", ph_units="eV")
#
#                plt.xlabel("energy (eV)")
#                plt.ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#                #plt.legend(["beam on cut, sprmc_th < " + str(self.sprmc_th),"on: {:.2f} eV, {:.0f}% tail".format(resoncut, tail_frac_oncut*100),
##                        "beam off","off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)],loc="best")
#                plt.title("{}, {} good chans\n{}".format(self.spill.runstr,self.data.n_good_channels(),linename))
#
#
#    # added by R.H. 20180620   
#    def plot_linefit_each(self,linename):
#        """
#        to plot each spectra 
#        """
#        ds0=self.data.first_good_dataset
#        elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
#        edges = np.arange(elo,ehi,1)
#        self.linefit_each_dict = {}
#        print "printing each spectra of %s starts ...."%(linename)
#        with PdfPages(self.get_output_base()+"run{:04d}_n{:04d}".format(self.first_pulse_runnum,self.noise_runnum)+ "_%s_nobeam.pdf"%linename) as pdf:
#            fig = plt.figure(figsize=(20,15))
#            coltmp=ds0.column_number
#            for ds in self.data:
#                col = ds.column_number
#                if col!=coltmp:
#                    fig.tight_layout()
#                    pdf.savefig()
#                    plt.close()
#                    fig = plt.figure(figsize=(20,15))
#                    coltmp=col
#                ax = fig.add_subplot(5,6,ds.row_number+1)
#                counts, centers = util.hist_ds(ds, edges, category={})
#                fitter = util.linefit(counts, edges, linename)
#                res = fitter.last_fit_params_dict["resolution"][0]
#                tail_frac = fitter.last_fit_params_dict["tail_frac"][0]
#                if res > 10: 
#                    ax.patch.set_facecolor('yellow')
#                    ax.patch.set_alpha(0.3)
#                fitter.plot(axis=ax,label=False, color="r", ph_units="eV")
#                ax.set_xlabel("energy (eV)")# + str(ds.channum))
#                ax.set_ylabel("counts per {:.1f} eV bin".format(centers[1]-centers[0]))
#                ax.legend(["no beam","no beam: {:.2f} eV, {:.0f}% tail".format(res, tail_frac*100)],loc="upper left",frameon=False,fontsize=10)
#                # ax.legend()
#                ax.set_title("chan{}".format(ds.channum))
##                print "{} :: chan{} -> no beam: {:.2f} eV, {:.0f}% tail".format(linename, ds.channum, res, tail_frac100), "| off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)
#                list = [res,tail_frac,0,0]
#                self.linefit_each_dict[ds.channum] = list
#
#            self.linefit_dict[linename] = self.linefit_each_dict
#            fig.tight_layout()
#            # fig.suptitle("run{:04d} %s (catecut%s)".format(self.first_pulse_runnum,linename,str(catecut)), fontsize=20)
#            # plt.subplots_adjust(top=0.90)
#            pdf.savefig()
#            plt.close()
#
#    # added by R.H. 20180620   
#    def plot_linefit_each_beam(self,linename):
#        """
#        to plot each spectra 
#        """
#        ds0=self.data.first_good_dataset
#        elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
#        edges = np.arange(elo,ehi,1)
#        self.linefit_each_dict = {}
#        print "printing each spectra of %s starts ...."%(linename)
#        with PdfPages(self.get_output_base()+"run{:04d}_n{:04d}".format(self.first_pulse_runnum,self.noise_runnum)+ "_%s.pdf"%linename) as pdf:
#            fig = plt.figure(figsize=(20,15))
#            coltmp=ds0.column_number
#            for ds in self.data:
#                col = ds.column_number
#                if col!=coltmp:
#                    fig.tight_layout()
#                    pdf.savefig()
#                    plt.close()
#                    fig = plt.figure(figsize=(20,15))
#                    coltmp=col
#                ax = fig.add_subplot(5,6,ds.row_number+1)
#                countson,  centers = util.hist_ds(ds, edges, category={"beam":"on"})
#                countsoff, centers = util.hist_ds(ds, edges, category={"beam":"off"})
#                fitteron = util.linefit(countson, edges, linename)
#                fitteroff = util.linefit(countsoff, edges, linename)
#                reson = fitteron.last_fit_params_dict["resolution"][0]
#                tail_frac_on = fitteron.last_fit_params_dict["tail_frac"][0]
#                resoff = fitteroff.last_fit_params_dict["resolution"][0]
#                tail_frac_off = fitteroff.last_fit_params_dict["tail_frac"][0]
#                if reson > 10: 
#                    ax.patch.set_facecolor('yellow')
#                    ax.patch.set_alpha(0.3)
#                fitteron.plot(axis=ax,label=False, color="r", ph_units="eV")
#                fitteroff.plot(axis=ax,label=False, color="b", ph_units="eV")
#                ax.set_xlabel("energy (eV) " + str(ds.channum))
#                ax.set_ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#                ax.legend(["beam on","on: {:.2f} eV, {:.0f}% tail".format(reson, tail_frac_on*100),"beam off","off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)],loc="upper left",frameon=False,fontsize=10)
#                # ax.legend()
#                ax.set_title("chan{}".format(ds.channum))
#                print "{} :: chan{} -> on: {:.2f} eV, {:.0f}% tail".format(linename, ds.channum, reson, tail_frac_on*100), "| off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)
#                list = [reson,tail_frac_on,resoff,tail_frac_off]
#                self.linefit_each_dict[ds.channum] = list
#
#            self.linefit_dict[linename] = self.linefit_each_dict
#            fig.tight_layout()
#            # fig.suptitle("run{:04d} %s (catecut%s)".format(self.first_pulse_runnum,linename,str(catecut)), fontsize=20)
#            # plt.subplots_adjust(top=0.90)
#            pdf.savefig()
#            plt.close()
#
#    def plot_linefit_beam_group_cut(self,linename):
#        elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
#        edges = np.arange(elo,ehi,1)
#        countson, centers =  util.hist_data(self.data, edges, category={"beam":"on","group_trigger":"far"})
#        countsoff, centers = util.hist_data(self.data, edges, category={"beam":"off","group_trigger":"far"})
#        fitteron = util.linefit(countson, edges, linename)
#        fitteroff = util.linefit(countsoff, edges, linename)
#        plt.figure()
#        fitteron.plot(axis=plt.gca(),label=False, color="r", ph_units="eV")
#        fitteroff.plot(axis=plt.gca(),label=False, color="b", ph_units="eV")
#        plt.xlabel("energy (eV)")
#        plt.ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#        reson = fitteron.last_fit_params_dict["resolution"][0]
#        tail_frac_on = fitteron.last_fit_params_dict["tail_frac"][0]
#        resoff = fitteroff.last_fit_params_dict["resolution"][0]
#        tail_frac_off = fitteroff.last_fit_params_dict["tail_frac"][0]
#        plt.legend(["beam on","on: {:.2f} eV, {:.0f}% tail".format(reson, tail_frac_on*100),"beam off","off: {:.2f} eV, {:.0f}% tail".format(resoff, tail_frac_off*100)],loc="best")
#        plt.title("{}, {} good chans\n{}".format(self.spill.runstr,self.data.n_good_channels(),linename))
#
#    def plot_overall_hist_beam(self):
#        self.plot_overall_hist(beam=True)
#
#    def plot_overall_hist(self,beam=False,elo=0,ehi=10000):
#        edges = np.arange(elo,ehi,4)
#        if beam:
#            countson, centers = util.hist_data(self.data, edges, category={"beam":"on"})
#            countsoff, centers = util.hist_data(self.data, edges, category={"beam":"off"})
#            plt.figure()
#            plt.plot(centers, countson, label="beam on")
#            plt.plot(centers, countsoff, label="beam off")
#        else:
#            counts, centers = util.hist_data(self.data, edges, category={})
#            plt.figure()
#            plt.plot(centers, counts, label="no beam")
#        plt.ylim(.9,plt.ylim()[1])
#        plt.yscale("log")
#        plt.xlabel("energy (eV)")
#        plt.ylabel("counts per {:.2f} eV bin".format(centers[1]-centers[0]))
#        plt.legend(loc="best")
#        plt.title("{}, {} good chans".format(self.spill.runstr,self.data.n_good_channels()))
#
#    def plot_each_resol_hist(self,linename,beam=False):
#        resol_file = self.get_output_base()+"run{:04d}_n{:04d}".format(self.first_pulse_runnum,self.noise_runnum)+ "_%s.csv"%(linename)
#        print "plot chan vs. resol by ", resol_file
#        if not os.path.exists(resol_file):
#            print "%s is not existed"%(resol_file)
#            print "please run 'def plot_each()' and 'def savetxt_each()' first"
#            sys.exit()
#
#        fin = open(resol_file,'r')
#        chan_list = []
#        resol_list1 = [] #on
#        resol_list2 = [] #off
#        for i, oneline in enumerate(fin):
#            if i == 0: continue
#            list = oneline.strip().split(',')
#            chan_list.append(int(list[0]))
#            resol_list1.append(float(list[1]))
#            if beam:                
#                resol_list2.append(float(list[3]))
#
#        fin.close()
#        self.chan_list = np.array(chan_list)
#        self.resol1 = np.array(resol_list1)
#        self.resol1_mean = np.mean(self.resol1[np.logical_not(np.isnan(self.resol1))])
#
#        gs = gridspec.GridSpec(1, 2, width_ratios=(3, 1))
#        ax = plt.subplot(gs[0, 0])
#        ax.set_title("run{:04d}, {} good chans\n {}".format(self.first_pulse_runnum,self.data.n_good_channels(),linename))
#
#        if not beam:
#            ax.plot(self.chan_list,self.resol1,color="b",marker=".",label="no beam")
#        else:
#            self.resol2 = np.array(resol_list2)
#            self.resol2_mean = np.mean(self.resol2[np.logical_not(np.isnan(self.resol1))])
#            ax.plot(self.chan_list,self.resol1,color="r",marker=".",label="beam on")
#            ax.plot(self.chan_list,self.resol2,color="b",marker=".",label="beam off")
#
#        ax.grid()
#        ax.legend()
#        ax.set_xlabel("chan")
#        ax.set_ylabel("FWHM energy resolution [eV]")
#        ax.set_xlim(0,480)
#        ax.set_ylim(0,20)
#
#        ax2 = plt.subplot(gs[0, 1],sharey=ax)
#        ax2.hist(self.resol1[np.logical_not(np.isnan(self.resol1))],bins=np.arange(3.,20.,0.25),histtype="step",color="r",orientation="horizontal")
#        if beam:
#            ax2.hist(self.resol2[np.logical_not(np.isnan(self.resol2))],bins=np.arange(3.,20.,0.25),histtype="step",color="b",orientation="horizontal")
#
#        ax2.tick_params(labelleft=False)
#        ax2.grid()
#        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.1)
#        if not beam:
#            plt.figtext(0.65,0.91,"resol mean =%.3f [eV]"%(self.resol1_mean))
#            plt.savefig(self.get_output_base()+"run{:04d}_n{:04d}".format(self.first_pulse_runnum,self.noise_runnum)+ "_resolhist_%s_nobeam.png"%linename)
#        else:
#            plt.figtext(0.65,0.95,"on_resol mean =%.3f [eV]"%(self.resol1_mean))
#            plt.figtext(0.65,0.91,"off_resol mean =%.3f [eV]"%(self.resol2_mean))
#            plt.savefig(self.get_output_base()+"run{:04d}_n{:04d}".format(self.first_pulse_runnum,self.noise_runnum)+ "_resolhist_%s.png"%linename)
#        plt.close()
#
#    def plot_each_resol_hist_beam(self,linename):
#        self.plot_each_resol_hist(linename,beam=True)
#
#    def plot_grptrig_pulse(self):
#        #grouptrig_csvfile = "../csv/grptrig_twocol.txt"
#        grouptrig_csvfile = self.GRTINFO
#        column_info_csvfile = self.COLUMN_INFO
#        print "start plot grptrig pulses..."        
#        (pulse_records_array, peak_region_max_all, peak_region_mean_all, 
#            peak_region_sum_all, energy_all, pretrig_mean_all, pretrig_rms_all, peak_value_all) = grp.get_grouptrig_pulse_info(self.data,self.pulse_cutarray, self.pulse_cut_category,self.GRTINFO,self.getMax,self.prim_chan_list)
#        
#        self.pulse_records_array = pulse_records_array
#        self.peak_region_max_all = peak_region_max_all
#        self.peak_region_mean_all = peak_region_mean_all
#        self.peak_region_sum_all = peak_region_sum_all
#        self.energy_all = energy_all
#        self.pretrig_mean_all = pretrig_mean_all
#        self.pretrig_rms_all = pretrig_rms_all
#        self.peak_value_all = peak_value_all
#
#        gdict = util.get_grtdict(grouptrig_csvfile)
#        col_dict = util.get_column_info(column_info_csvfile)
#        # pprint.pprint(pulse_records_array)
#        for i, chan in enumerate(self.prim_chan_list):
#        # for i, chan in enumerate(self.prim_chan_list[0]):
#            grtchs = gdict.get(chan)
#            grtchs.insert(0,chan)
#            grtchs = [x for x in grtchs if not x in self.badchans]
#            if self.energy_cut == True and self.cut_category != "None":
#                figname = self.get_output_base()+"run{:04d}".format(self.first_pulse_runnum)+ "_grppulse_%s_Ecut_%s_%s_%s_to_%s.pdf"%(str(chan),str(self.energy_cut),self.pulse_cut_category, str(self.cutMin), str(self.cutMax))
#                pagetitle = "run%04d Chan%d w/ %s cut (%s - %s) energy_cut = %s (%.2f - %.2f)" %(self.first_pulse_runnum, chan,self.pulse_cut_category, str(self.cutMin), str(self.cutMax), 
#                                                  str(self.energy_cut), float(self.energyMin), float(self.energyMax))                                     
#            elif self.energy_cut == True and self.cut_category == "None":
#                figname = self.get_output_base()+"run{:04d}".format(self.first_pulse_runnum)+ "_grppulse_%s_Ecut_%s_to_%s.pdf"%(str(chan),str(self.energy_cut), str(self.energyMin), str(self.energyMax))
#                pagetitle = "run%04d Chan%d w/o any cuts energy_cut = %s (%.2f - %.2f)"%(self.first_pulse_runnum,chan,str(self.energy_cut), float(self.energyMin), float(self.energyMax))       
#            elif self.energy_cut == False and self.cut_category != "None":
#                figname = self.get_output_base()+"run{:04d}".format(self.first_pulse_runnum)+ "_grppulse_%s_Ecut_%s_%s_%s_to_%s.pdf"%(str(chan),str(self.energy_cut),self.pulse_cut_category, str(self.cutMin), str(self.cutMax))
#                pagetitle = "run%04d Chan%d w/ %s cut (%s - %s) energy_cut = %s "%(self.first_pulse_runnum, chan,self.pulse_cut_category, str(self.cutMin), str(self.cutMax), str(self.energy_cut)) 
#            elif self.energy_cut == False and self.cut_category == "None":
#                figname = self.get_output_base()+"run{:04d}".format(self.first_pulse_runnum)+ "_grppulse_%s.pdf"%(str(chan))
#                pagetitle = "run%04d Chan%d w/o any cuts"%(self.first_pulse_runnum, chan)
#
#            self.figname = figname
#
#            with PdfPages(figname) as pdf:
#                fig = plt.figure(figsize=(20,15))
#                if not self.grouptrigmax: ev = 0
#                for j in xrange(len(pulse_records_array[i])):
#                    print "Dataset No.%d, plot %d pulses"%(j+1,len(pulse_records_array[i][j]))
#                    
#                    if self.grouptrigmax: ev = 0
#
#                    for one_grpchan in grtchs:
#                        # ev += 1
#                        if ev%20 == 0 and ev == 0:
#                            k = 1
#                        else:
#                            k = ev%20 +1
#
#                        if self.grouptrigmax:
#                        #if not self.calibration_runnum is None:
#                            # if k%5 ==  1 and ev==0: color = "b"
#                            # elif k%5 == 1 and ev != 0: color = "k"
#                            # elif k%5 == 2: color = "r"
#                            # elif k%5 == 3: color = "g"
#                            # elif k%5 == 4: color = "c"
#                            # elif k%5 == 0: color = "m"
#                            col_num = tesmap.getCol(one_grpchan)
#                            if col_num == 0: col_name = "AX"; color = "k"
#                            elif col_num == 1: col_name = "AY"; color = "r"
#                            elif col_num == 2: col_name = "BX"; color = "g"
#                            elif col_num == 3: col_name = "BY"; color = "c"
#                            elif col_num == 4: col_name = "CX"; color = "m"
#                            elif col_num == 5: col_name = "CY"; color = "y"
#                            elif col_num == 6: col_name = "DX"; color = "orange"
#                            elif col_num == 7: col_name = "DY"; color = "hotpink"
#
#                            if one_grpchan == chan: color = "b"
#                        else:
#                            if ev%len(grtchs) == 0: color = "b"
#                            elif ev%len(grtchs) == 1: color = "r"
#                            elif ev%len(grtchs) == 2 : color = "g"
#                            elif ev%len(grtchs) == 3 : color = "c"
#                            elif ev%len(grtchs) == 4 : color = "m"
#
#
#                        # if not self.calibration_runnum is None:
#                        #     color = cm.jet(k/20)
#
#                        distance = tesmap.getDistance(chan,one_grpchan)
#
#                        ax = fig.add_subplot(4,5,k)
#                        if one_grpchan in self.badchans: continue
#                        ax.plot(pulse_records_array[i][j][one_grpchan],color=color, 
#                            label="energy=%.2f\npr_max=%.2f\npr_mean=%.2f\npr_sum=%.2f\npretrig_mean=%.2f\npeak_value=%.2f"
#                            %(energy_all[i][j][one_grpchan],peak_region_max_all[i][j][one_grpchan], 
#                              peak_region_mean_all[i][j][one_grpchan], peak_region_sum_all[i][j][one_grpchan], 
#                              pretrig_mean_all[i][j][one_grpchan], peak_value_all[i][j][one_grpchan]))
#                        ax.set_title("chan%d(d=%.2fum, RS=%s)"%(one_grpchan, distance, col_dict[str(one_grpchan)]))
#                        ax.legend(loc="upper right",fontsize=8)
#                        if self.cut_category == "sec_pr_mean" and one_grpchan != chan and self.grouptrigmax != None:
#                            if self.cutMin == None:
#                                pr_mean_value = np.array(peak_region_mean_all[i][j].values())
#                                # pr_mean_value = pr_mean_value[np.where(pr_mean_value < self.cutMax)]
#                                nearest_ind = (np.abs(pr_mean_value - self.cutMax)).argmin()
#                                # print peak_region_mean_all[i][j].values()[nearest_ind]
#                            
#                            if self.cutMax == None:
#                                pr_mean_value = np.array(peak_region_mean_all[i][j].values())
#                                # pr_mean_value = pr_mean_value[np.where(pr_mean_value > self.cutMin)]
#                                nearest_ind = (np.abs(pr_mean_value - self.cutMin)).argmin()
#                                
#                            if peak_region_mean_all[i][j][one_grpchan] == peak_region_mean_all[i][j].values()[nearest_ind]: 
#                                ax.patch.set_facecolor('yellow')
#                                ax.patch.set_alpha(0.3)
#                            # elif peak_region_mean_all[i][j][one_grpchan] > self.cutMin or peak_region_mean_all[i][j][one_grpchan] < self.cutMax:
#                            #     ax.patch.set_facecolor('yellow')
#                            #     ax.patch.set_alpha(0.3)
#                        elif self.grouptrigmax and ev == 0:
#                            ax.patch.set_facecolor('cyan')
#                            ax.patch.set_alpha(0.3)
#                        elif len(grtchs) < 5 and one_grpchan == chan:
#                            ax.patch.set_facecolor('cyan')
#                            ax.patch.set_alpha(0.3)
#                                
#                        ev += 1
#                        fig.tight_layout()
#
#                        if self.grouptrigmax:
#                        #if not self.calibration_runnum is None:
#                            if k == 20 and ev != 0 or ev == len(pulse_records_array[i][j]):
#                                if self.energy_cut:
#                                    fig.suptitle(pagetitle,fontsize=20)
#                                else:
#                                    fig.suptitle(pagetitle,fontsize=20)
#                                plt.subplots_adjust(top=0.93)
#                                #plt.savefig(self.get_output_base()+"run{:04d}_n{:04d}".format(self.first_pulse_runnum,self.noise_runnum)+ "_grppulse_%s_%d.png"%(str(chan),ev))
#                                pdf.savefig()
#                                plt.clf()
#                    
#                    if not self.grouptrigmax:
#                    #if self.calibration_runnum is None:
#                        if k == 20 and ev != 0:
#                            fig.suptitle(pagetitle,fontsize=20)
#                            plt.subplots_adjust(top=0.93)
#                            pdf.savefig()
#                            plt.clf()
#
#                if self.grouptrigmax:
#                # if not self.calibration_runnum is None:
#                    #plt.figure(figsize=(20,15))
#                    ax2 = tesmap.getPlotRectPlane(1,1,1)
#                    for ii in xrange(1,481,2):
#                        col_num = tesmap.getCol(ii)
#                        if col_num == 0: col_name = "AX"; fcolor = "k"
#                        elif col_num == 1: col_name = "AY"; fcolor = "r"
#                        elif col_num == 2: col_name = "BX"; fcolor = "g"
#                        elif col_num == 3: col_name = "BY"; fcolor = "c"
#                        elif col_num == 4: col_name = "CX"; fcolor = "m"
#                        elif col_num == 5: col_name = "CY"; fcolor = "y"
#                        elif col_num == 6: col_name = "DX"; fcolor = "orange"
#                        elif col_num == 7: col_name = "DY"; fcolor = "hotpink"
#                        if ii == chan:
#                            chanpos = tesmap.getRectTESPos(ii,facecolor='b',edgecolor='red',alpha=0.5)
#                            ax2.add_patch(chanpos)
#                        else:
#                            chanpos1 = tesmap.getRectTESPos(ii,facecolor=fcolor,edgecolor=fcolor,alpha=0.4)
#                            ax2.add_patch(chanpos1)
#                        tesmap.writeTextTESPos(ax2,ii,ii,fontsize=15)
#                    ax2.set_title("pulse color and TES position for groupTrig MAX")
#                    ax2.set_xlabel("Position X [$\mu$m]")
#                    ax2.set_ylabel("Position Y [$\mu$m]")
#                    plt.tight_layout()
#                    pdf.savefig()
#                    plt.clf()
#            
#            print "fig saved at ",figname
#            plt.close()
#
#    def plot_grptrig_peak_value_map(self):
#        #grouptrig_csvfile = "../csv/grptrig_twocol.txt"
#        grouptrig_csvfile = self.GRTINFO
#        column_info_csvfile = self.COLUMN_INFO
#        print "start plot grptrig peak_value map..."
#        # (pulse_records_array, peak_region_max_all, peak_region_mean_all, 
#        #     peak_region_sum_all, energy_all, pretrig_mean_all, pretrig_rms_all, peak_value_all) = khe.get_grouptrig_pulse_info(self.data,self.pulse_cutarray, self.pulse_cut_category,self.GRTINFO,self.getMax,self.prim_chan_list[0])
#        peak_value_all = self.peak_value_all
#        gdict = util.get_grtdict(grouptrig_csvfile)
#        col_dict = util.get_column_info(column_info_csvfile)
#        # pprint.pprint(pulse_records_array)
#        for i, chan in enumerate(self.prim_chan_list):
#        # for i, chan in enumerate(self.prim_chan_list[0]):
#            grtchs = gdict.get(chan)
#            grtchs.insert(0,chan)
#            grtchs = [x for x in grtchs if not x in self.badchans]
#            #figname = self.get_output_base()+"run{:04d}".format(self.first_pulse_runnum)+ "_grpPVmap_%s_Ecut_%s_%s_%s_to_%s.pdf"%(str(chan),str(self.energy_cut),self.cut_category, str(self.cutMin), str(self.cutMax))
#            figname = self.figname.replace("pulse","PVmap")
#            with PdfPages(figname) as pdf:
#                plt.figure(figsize=(8.27, 8.27))
#                #ax = tesmap.getPlotRectPlane(1,1,1,)
#                ev = 0
#                for j in xrange(len(peak_value_all[i])):
#                    ax = tesmap.getPlotRectPlane(1,1,1)
#                    for ii in xrange(1,481,2):
#                        tesmap.writeTextTESPos(ax,ii,ii,fontsize=8)
#                    # print "len(pulse_records_array[i][j])=",len(pulse_records_array[i][j])
#                    pv_array = np.array(peak_value_all[i][j].values())
#                    ax = tesmap.getPlotDensityPlane(1,1,1,chans=grtchs[1::],resols=pv_array[np.where(pv_array!=pv_array.max())])
#                    for jj,one_grpchan in enumerate(grtchs):
#                        tesmap.writeFloatTESPos(ax,one_grpchan,float(peak_value_all[i][j][one_grpchan]), fontsize=6, yoffset=250.0, color='red')
#                        if one_grpchan == chan:
#                            chanpos1 = tesmap.getRectTESPos(one_grpchan,facecolor='red',edgecolor='red',alpha=0.5)
#                            ax.add_patch(chanpos1)
#                    plt.xlabel("Position x [$\mu$m]")
#                    plt.ylabel("Position y [$\mu$m]")
#                    plt.title("run{:04d} peak value of each pixels\nprimary pulse peak value = {:.1f}".format(self.first_pulse_runnum,float(peak_value_all[i][j][chan])))
#                    pdf.savefig()
#                    plt.clf()
#            print "fig saved at ",figname
#            plt.close()
