import os
import pylab as plt
import pandas as pd
import numpy as np
import csv
from matplotlib.backends.backend_pdf import PdfPages

outdir="./output/"

#runs=[36,38,40,43,44,54,61,63,66,67,76,78,84,85,87,88,97,98,100]# all
runs=[36,38,43,54,61,67,76,84,88,98,100]# MnKa

files=[outdir+"run%d_resol.csv"%run for run in runs]
chs=[]
resols=[]
resols_err=[]
for run, f in zip(runs,files):
    if not os.path.isfile(f):
        print "%s does not exist"%f
        continue
    df = pd.read_csv(f,header=None)
    chs.append(df.iloc[:,0].tolist())
    resols.append(df.iloc[:,1].tolist())
    resols_err.append(df.iloc[:,2].tolist())

    
plt.close('all')
plt.ion()
  

fwhmth=10.
#howmanyover=0
#howmanyover=1
#howmanyover=2
howmanyover=3

ncols=8
nrows=30
cols=["AX","AY","BX","BY","CX","CY","DX","DY"]

nch=0
goodch=[]
LAST=False
foutname="resols_MnKa_th%d_n%d"%(int(fwhmth),howmanyover)
fout = open(outdir+foutname+".csv", 'w')
writer = csv.writer(fout, lineterminator='\n')
pdfname=outdir+foutname+".pdf"
with PdfPages(pdfname) as pdf:
    for j in xrange(ncols):
        fig = plt.figure()
        fig.set_size_inches(11.69,8.27)
        fig.suptitle(foutname+"_%s"%cols[j])
        ichs = [j*nrows*2+ich*2-1 for ich in xrange(1,nrows+1,1)]
        for k in xrange(6):
            if LAST: break
            ax = plt.subplot(3,2,k+1)
            #for ch in chs[0][5*k+nrows*j:5*(k+1)+nrows*j]:
            for ch in ichs[5*k:5*(k+1)]:
                if not ch in chs[0]: continue
                if ch==479: LAST=True
                nch+=1
                x=runs
                y=np.zeros(len(runs),dtype=np.float32)
                yerr=np.zeros(len(runs),dtype=np.float32)
                for i,run in enumerate(runs):
                    y[i]=resols[i][chs[i].index(ch)]
                    yerr[i]=resols_err[i][chs[i].index(ch)]
                    if y[i]>100 or yerr[i]>100:
                        y[i]=50.
                        yerr[i]=0.
                label='ch%d bad'%ch
                if len(np.where(y>fwhmth)[0])<=howmanyover:
                    label='ch%d good'%ch
                    goodch.append(["%d"%ch])
                    print ch
                ax.errorbar(x,y,yerr=yerr,fmt='o',elinewidth=1,label=label)
            ax.set_ylim(1,20)
            ax.set_xlabel("run number")
            ax.set_ylabel("FWHM eV at MnKAlpha")
            ax.legend(loc="upper right")
        fig.tight_layout()
        fig.subplots_adjust(top=0.92)
        pdf.savefig()
        plt.close()
        

print nch,len(goodch),float(len(goodch))/nch
writer.writerows(goodch)
fout.close()

