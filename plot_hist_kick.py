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
cut=""
fnames = [ROOTDIR+"/run0098/run0098_noi0086_mass_2019%s.root"%cut,
          ROOTDIR+"/run0100/run0100_noi0086_mass_2019%s.root"%cut]
tags=["kickoff","kickon"]

# good channel selection
fwhmth=10.
howmany_n=3
#howmany_n=1
fgoodch=outdir+"resols_MnKa_th%d_n%d.csv"%(int(fwhmth),howmany_n)
df = pd.read_csv(fgoodch,header=None)
goodchs = df.iloc[:,0].tolist()
print "good channels: ", goodchs

fs = [ROOT.TFile.Open(fname,"read") for fname in fnames]
foutname=outdir+"fe55_kicker_onoff_n%d.root"%howmany_n
fout = ROOT.TFile.Open(foutname,"recreate")

for f,tag in zip(fs,tags):
    f.cd()
    t = f.Get("chanall")

    fout.cd()
    hene_off   = ROOT.TH1F("hene_off_%s"%tag,"hene_off_%s"%tag,15000,0,15000)
    hene_offFill=hene_off.Fill
    set_hist(hene_off,xtitle="Energy [eV]",
             ytitle="Counts / %.1f eV"%(hene_off.GetBinWidth(1)))

    # ---- loop ----
    chtmp=0
    for e in t:
        if not e.ch in goodchs: continue
        if e.good!=1: continue
        elif e.primary!=1: continue
        if chtmp<e.ch:
            chtmp=e.ch
            print tag, chtmp
        hene_offFill(e.energy)
        
    hene_off.Write()  
    hene_off.Delete()  
    f.Close()

fout.Close()
# ---------------
