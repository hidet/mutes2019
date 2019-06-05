import os
import pylab as plt
import pandas as pd
import numpy as np
import csv
import ROOT


def event_check(e):
    if e.good!=1: return False
    elif e.primary!=1: return False
    elif e.jbrsc!=1: return False
    elif e.sprmc==1: return True
    else: return False


def set_hist(h,xtitle="",ytitle=""):
    h.Sumw2()
    h.SetTitle(h.GetName())
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)


outdir="./output/"
ROOTDIR="%s/dumproot"%(os.environ['MUTESDATADIR'])
cut="_spilloff_sprmcon_jbrscon_prime"
fnames = [ROOTDIR+"/run0036/run0036_noi0035_mass_2019%s.root"%cut,
          ROOTDIR+"/run0058/run0058_noi0055_mass_2019_58_59%s.root"%cut,
          ROOTDIR+"/run0079/run0079_noi0077_mass_2019_79_80_81%s.root"%cut]
tags=["Ne09atm","Ne04atm","Ne01atm"]

# good channel selection
fwhmth=10.
howmany_n=3
#howmany_n=1
fgoodch=outdir+"resols_MnKa_th%d_n%d.csv"%(int(fwhmth),howmany_n)
df = pd.read_csv(fgoodch,header=None)
goodchs = df.iloc[:,0].tolist()
print "good channels: ", goodchs

fs = [ROOT.TFile.Open(fname,"read") for fname in fnames]
foutname=outdir+"MuNeXray_n%d.root"%howmany_n
fout = ROOT.TFile.Open(foutname,"recreate")

for f,tag in zip(fs,tags):
    f.cd()
    t = f.Get("chanall")

    fout.cd()
    hene_off   = ROOT.TH1F("hene_off_%s"%tag,"hene_off_%s"%tag,15000,0,15000)
    hene_bon   = ROOT.TH1F("hene_bon_%s"%tag,"hene_bon_%s"%tag,15000,0,15000)
    hene_ton   = ROOT.TH1F("hene_ton_%s"%tag,"hene_ton_%s"%tag,15000,0,15000)
    htime_all  = ROOT.TH1F("htime_all_%s"%tag,"htime_all_%s"%tag,600,0,150)
    htime_mune = ROOT.TH1F("htime_mune_%s"%tag,"htime_mune_%s"%tag,600,0,150)
    h2         = ROOT.TH2F("h2_%s"%tag,"h2_%s"%tag,600,0,150,450,5800,6700)

    hene_offFill=hene_off.Fill
    hene_bonFill=hene_bon.Fill
    hene_tonFill=hene_ton.Fill
    htime_allFill=htime_all.Fill
    htime_muneFill=htime_mune.Fill
    h2Fill=h2.Fill

    set_hist(hene_off,xtitle="Energy [eV]",
             ytitle="Counts / %.1f eV"%(hene_off.GetBinWidth(1)))
    set_hist(hene_bon,xtitle="Energy [eV]",
             ytitle="Counts / %.1f eV"%(hene_bon.GetBinWidth(1)))
    set_hist(hene_ton,xtitle="Energy [eV]",
             ytitle="Counts / %.1f eV"%(hene_ton.GetBinWidth(1)))
    set_hist(htime_all,xtitle="Time (1ch=240ns)",
             ytitle="Counts / %.1f ch"%(htime_all.GetBinWidth(1)))
    set_hist(htime_mune,xtitle="Time (1ch=240ns)",
             ytitle="Counts / %.1f ch"%(htime_mune.GetBinWidth(1)))
    set_hist(h2,xtitle="Time (1ch=240ns)",ytitle="Energy [eV]")


    # ---- loop ----
    chtmp=0
    for e in t:
        if not e.ch in goodchs: continue
        if not event_check(e): continue
        if chtmp<e.ch:
            chtmp=e.ch
            print tag, chtmp
        if e.beam==0:# beam off
            hene_offFill(e.energy)
        elif e.beam==1:# beam on
            hene_bonFill(e.energy)
            htime_allFill(e.row_next_extrig_nrp)
            if e.row_next_extrig_nrp>=0 and e.row_next_extrig_nrp<150:
                h2Fill(e.row_next_extrig_nrp,e.energy)
            if e.energy>=6260 and e.energy<6340:
                htime_muneFill(e.row_next_extrig_nrp)
            #if e.row_next_extrig_nrp>=65 and e.row_next_extrig_nrp<75:
            if e.row_next_extrig_nrp>=64 and e.row_next_extrig_nrp<74:
                hene_tonFill(e.energy)
            
    hene_off.Write()  
    hene_bon.Write()    
    hene_ton.Write()    
    htime_all.Write()   
    htime_mune.Write()  
    h2.Write()                      

    hene_off.Delete()  
    hene_bon.Delete()    
    hene_ton.Delete()    
    htime_all.Delete()   
    htime_mune.Delete()  
    h2.Delete()

    f.Close()

fout.Close()
# ---------------
