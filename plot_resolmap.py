import os
import pylab as plt
import pandas as pd
import numpy as np
import csv
import tesmap_forgrptrig as tesmap
from matplotlib.backends.backend_pdf import PdfPages

outdir="./output/"

#runs=[36,38,40,43,44,54,61,63,66,67,76,78,84,85,87,88,97,98,100]# all
runs=[36,38,43,54,61,67,76,84,88,98,100]# MnKa

files=[outdir+"TMU_2019G/run%d_resol.csv"%run for run in runs]
    
plt.close('all')
plt.ion()

fwhmth=10.
#howmanyover=0
#howmanyover=1
#howmanyover=2
howmanyover=3

foutname="resolsmap_MnKa_th%d_n%d"%(int(fwhmth),howmanyover)
pdfname=outdir+foutname+".pdf"
with PdfPages(pdfname) as pdf:
    for f in files:
        print f
        tesmap.plotDensityPlaneFromCSV(f)
        pdf.savefig()
        plt.close()

