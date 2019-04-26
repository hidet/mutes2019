import numpy as np
from math import pow, sqrt
#import matplotlib
#import matplotlib.pyplot as plt
import pylab as plt
import matplotlib.cm as cm
import matplotlib as mpl

spacing = 520

position=[
    [  1,      290   ,   -3470,    0 ],
	[  3,      290   ,   -4510,    1 ],
	[  5,      290   ,   -1910,    2 ],
	[  7,      290   ,   -870 ,    3 ],
	[  9,      870   ,   -2430,    4 ],
	[ 11,      870   ,   -3470,    5 ],
	[ 13,      870   ,   -4510,    6 ],
	[ 15,      870   ,   -870 ,    7 ],
	[ 17,      1450  ,   -2030,    8 ],
	[ 19,      1450  ,   -2550,    9 ],
	[ 21,      1450  ,   -3590,    10],
	[ 23,      3190  ,   -3190,    11],
	[ 25,      2610  ,   -3710,    12],
	[ 27,      2030  ,   -3590,    13],
	[ 29,      2030  ,   -2550,    14],
	[ 31,      290   ,   -2950,    15],
	[ 33,      290   ,   -3990,    16],
	[ 35,      290   ,   -2430,    17],
	[ 37,      290   ,   -1390,    18],
	[ 39,      870   ,   -3990,    19],
	[ 41,      870   ,   -1910,    20],
	[ 43,      870   ,   -1390,    21],
	[ 45,      2030  ,   -2030,    22],
	[ 47,      1450  ,   -3070,    23],
	[ 49,      1450  ,   -4110,    24],
	[ 51,      870   ,   -2950,    25],
	[ 53,      3190  ,   -3710,    26],
	[ 55,      2610  ,   -3190,    27],
	[ 57,      2030  ,   -3070,    28],
	[ 59,      2030  ,   -4110,    29],
	[ 61,      -2030 ,   -3650,    0 ],
	[ 63,      -2030 ,   -2610,    1 ],
	[ 65,      -2610 ,   -3130,    2 ],
	[ 67,      -3190 ,   -3710,    3 ],
	[ 69,      -870  ,   -1970,    4 ],
	[ 71,      -1450 ,   -4050,    5 ],
	[ 73,      -1450 ,   -3010,    6 ],
	[ 75,      -1450 ,   -1450,    7 ],
	[ 77,      -870  ,   -4570,    8 ],
	[ 79,      -870  ,   -4050,    9 ],
	[ 81,      -870  ,   -3010,    10],
	[ 83,      -290  ,   -810 ,    11],
	[ 85,      -290  ,   -3930,    12],
	[ 87,      -290  ,   -2890,    13],
	[ 89,      -290  ,   -1850,    14],
	[ 91,      -2030 ,   -4170,    15],
	[ 93,      -2030 ,   -3130,    16],
	[ 95,      -2610 ,   -2610,    17],
	[ 97,      -2610 ,   -3650,    18],
	[ 99,      -1450 ,   -3530,    19],
	[101,      -1450 ,   -2490,    20],
	[103,      -1450 ,   -1970,    21],
	[105,      -290  ,   -1330,    22],
	[107,      -870  ,   -3530,    23],
	[109,      -870  ,   -2490,    24],
	[111,      -870  ,   -1450,    25],
	[113,      -290  ,   -290 ,    26],
	[115,      -290  ,   -4450,    27],
	[117,      -290  ,   -3410,    28],
	[119,      -290  ,   -2370,    29],
	[121,      -3470 ,   -290 ,    0 ],
	[123,      -4510 ,   -290 ,    1 ],
	[125,      -1910 ,   -290 ,    2 ],
	[127,      -870  ,   -290 ,    3 ],
	[129,      -2430 ,   -870 ,    4 ],
	[131,      -3470 ,   -870 ,    5 ],
	[133,      -4510 ,   -870 ,    6 ],
	[135,      -870  ,   -870 ,    7 ],
	[137,      -2030 ,   -1450,    8 ],
	[139,      -2550 ,   -1450,    9 ],
	[141,      -3590 ,   -1450,    10],
	[143,      -3190 ,   -3190,    11],
	[145,      -3710 ,   -2610,    12],
	[147,      -3590 ,   -2030,    13],
	[149,      -2550 ,   -2030,    14],
	[151,      -2950 ,   -290 ,    15],
	[153,      -3990 ,   -290 ,    16],
	[155,      -2430 ,   -290 ,    17],
	[157,      -1390 ,   -290 ,    18],
	[159,      -3990 ,   -870 ,    19],
	[161,      -1910 ,   -870 ,    20],
	[163,      -1390 ,   -870 ,    21],
	[165,      -2030 ,   -2030,    22],
	[167,      -3070 ,   -1450,    23],
	[169,      -4110 ,   -1450,    24],
	[171,      -2950 ,   -870 ,    25],
	[173,      -3710 ,   -3190,    26],
	[175,      -3190 ,   -2610,    27],
	[177,      -3070 ,   -2030,    28],
	[179,      -4110 ,   -2030,    29],
	[181,      -3650 ,   2030 ,    0 ],
	[183,      -2610 ,   2030 ,    1 ],
	[185,      -3130 ,   2610 ,    2 ],
	[187,      -3710 ,   3190 ,    3 ],
	[189,      -1970 ,   870  ,    4 ],
	[191,      -4050 ,   1450 ,    5 ],
	[193,      -3010 ,   1450 ,    6 ],
	[195,      -1450 ,   1450 ,    7 ],
	[197,      -4570 ,   870  ,    8 ],
	[199,      -4050 ,   870  ,    9 ],
	[201,      -3010 ,   870  ,    10],
	[203,      -810  ,   290  ,    11],
	[205,      -3930 ,   290  ,    12],
	[207,      -2890 ,   290  ,    13],
	[209,      -1850 ,   290  ,    14],
	[211,      -4170 ,   2030 ,    15],
	[213,      -3130 ,   2030 ,    16],
	[215,      -2610 ,   2610 ,    17],
	[217,      -3650 ,   2610 ,    18],
	[219,      -3530 ,   1450 ,    19],
	[221,      -2490 ,   1450 ,    20],
	[223,      -1970 ,   1450 ,    21],
	[225,      -1330 ,   290  ,    22],
	[227,      -3530 ,   870  ,    23],
	[229,      -2490 ,   870  ,    24],
	[231,      -1450 ,   870  ,    25],
	[233,      -290  ,   290  ,    26],
	[235,      -4450 ,   290  ,    27],
	[237,      -3410 ,   290  ,    28],
	[239,      -2370 ,   290  ,    29],
	[241,      -290  ,   3470 ,    0 ],
	[243,      -290  ,   4510 ,    1 ],
	[245,      -290  ,   1910 ,    2 ],
	[247,      -290  ,   870  ,    3 ],
	[249,      -870  ,   2430 ,    4 ],
	[251,      -870  ,   3470 ,    5 ],
	[253,      -870  ,   4510 ,    6 ],
	[255,      -870  ,   870  ,    7 ],
	[257,      -1450 ,   2030 ,    8 ],
	[259,      -1450 ,   2550 ,    9 ],
	[261,      -1450 ,   3590 ,    10],
	[263,      -3190 ,   3190 ,    11],
	[265,      -2610 ,   3710 ,    12],
	[267,      -2030 ,   3590 ,    13],
	[269,      -2030 ,   2550 ,    14],
	[271,      -290  ,   2950 ,    15],
	[273,      -290  ,   3990 ,    16],
	[275,      -290  ,   2430 ,    17],
	[277,      -290  ,   1390 ,    18],
	[279,      -870  ,   3990 ,    19],
	[281,      -870  ,   1910 ,    20],
	[283,      -870  ,   1390 ,    21],
	[285,      -2030 ,   2030 ,    22],
	[287,      -1450 ,   3070 ,    23],
	[289,      -1450 ,   4110 ,    24],
	[291,      -870  ,   2950 ,    25],
	[293,      -3190 ,   3710 ,    26],
	[295,      -2610 ,   3190 ,    27],
	[297,      -2030 ,   3070 ,    28],
	[299,      -2030 ,   4110 ,    29],
	[301,      2030  ,   3650 ,    0 ],
	[303,      2030  ,   2610 ,    1 ],
	[305,      2610  ,   3130 ,    2 ],
	[307,      3190  ,   3710 ,    3 ],
	[309,      870   ,   1970 ,    4 ],
	[311,      1450  ,   4050 ,    5 ],
	[313,      1450  ,   3010 ,    6 ],
	[315,      1450  ,   1450 ,    7 ],
	[317,      870   ,   4570 ,    8 ],
	[319,      870   ,   4050 ,    9 ],
	[321,      870   ,   3010 ,    10],
	[323,      290   ,   810  ,    11],
	[325,      290   ,   3930 ,    12],
	[327,      290   ,   2890 ,    13],
	[329,      290   ,   1850 ,    14],
	[331,      2030  ,   4170 ,    15],
	[333,      2030  ,   3130 ,    16],
	[335,      2610  ,   2610 ,    17],
	[337,      2610  ,   3650 ,    18],
	[339,      1450  ,   3530 ,    19],
	[341,      1450  ,   2490 ,    20],
	[343,      1450  ,   1970 ,    21],
	[345,      290   ,   1330 ,    22],
	[347,      870   ,   3530 ,    23],
	[349,      870   ,   2490 ,    24],
	[351,      870   ,   1450 ,    25],
	[353,      290   ,   290  ,    26],
	[355,      290   ,   4450 ,    27],
	[357,      290   ,   3410 ,    28],
	[359,      290   ,   2370 ,    29],
	[361,      3470  ,   290  ,    0 ],
	[363,      4510  ,   290  ,    1 ],
	[365,      1910  ,   290  ,    2 ],
	[367,      870   ,   290  ,    3 ],
	[369,      2430  ,   870  ,    4 ],
	[371,      3470  ,   870  ,    5 ],
	[373,      4510  ,   870  ,    6 ],
	[375,      870   ,   870  ,    7 ],
	[377,      2030  ,   1450 ,    8 ],
	[379,      2550  ,   1450 ,    9 ],
	[381,      3590  ,   1450 ,    10],
	[383,      3190  ,   3190 ,    11],
	[385,      3710  ,   2610 ,    12],
	[387,      3590  ,   2030 ,    13],
	[389,      2550  ,   2030 ,    14],
	[391,      2950  ,   290  ,    15],
	[393,      3990  ,   290  ,    16],
	[395,      2430  ,   290  ,    17],
	[397,      1390  ,   290  ,    18],
	[399,      3990  ,   870  ,    19],
	[401,      1910  ,   870  ,    20],
	[403,      1390  ,   870  ,    21],
	[405,      2030  ,   2030 ,    22],
	[407,      3070  ,   1450 ,    23],
	[409,      4110  ,   1450 ,    24],
	[411,      2950  ,   870  ,    25],
	[413,      3710  ,   3190 ,    26],
	[415,      3190  ,   2610 ,    27],
	[417,      3070  ,   2030 ,    28],
	[419,      4110  ,   2030 ,    29],
	[421,      3650  ,   -2030,    0 ],
	[423,      2610  ,   -2030,    1 ],
	[425,      3130  ,   -2610,    2 ],
	[427,      3710  ,   -3190,    3 ],
	[429,      1970  ,   -870 ,    4 ],
	[431,      4050  ,   -1450,    5 ],
	[433,      3010  ,   -1450,    6 ],
	[435,      1450  ,   -1450,    7 ],
	[437,      4570  ,   -870 ,    8 ],
	[439,      4050  ,   -870 ,    9 ],
	[441,      3010  ,   -870 ,    10],
	[443,      810   ,   -290 ,    11],
	[445,      3930  ,   -290 ,    12],
	[447,      2890  ,   -290 ,    13],
	[449,      1850  ,   -290 ,    14],
	[451,      4170  ,   -2030,    15],
	[453,      3130  ,   -2030,    16],
	[455,      2610  ,   -2610,    17],
	[457,      3650  ,   -2610,    18],
	[459,      3530  ,   -1450,    19],
	[461,      2490  ,   -1450,    20],
	[463,      1970  ,   -1450,    21],
	[465,      1330  ,   -290 ,    22],
	[467,      3530  ,   -870 ,    23],
	[469,      2490  ,   -870 ,    24],
	[471,      1450  ,   -870 ,    25],
	[473,      290   ,   -290 ,    26],
	[475,      4450  ,   -290 ,    27],
	[477,      3410  ,   -290 ,    28],
	[479,      2370  ,   -290 ,    29]]


pos = np.array(position)


twradd = [7,5,3,1,15,13,11,9,
          23,19,17,31,29,27,25,
          6,4,2,0,14,12,8,
          22,20,18,16,30,28,26,24]

pixadd = [6,4,2,0,16,14,12,8,
          22,20,18,30,28,26,24,
          7,5,3,1,13,11,9,
          23,19,17,15,31,29,27,25]



def getCol(ch):
    if   ch>=1   and ch<=59: return 0 
    elif ch>=61  and ch<=119:return 1
    elif ch>=121 and ch<=179:return 2
    elif ch>=181 and ch<=239:return 3
    elif ch>=241 and ch<=299:return 4
    elif ch>=301 and ch<=359:return 5
    elif ch>=361 and ch<=419:return 6
    elif ch>=421 and ch<=479:return 7
    else: return -1


# add for MUTES
def getCollist(ch):
    if   ch>=1   and ch<=59: return np.arange(1,   59 + 1 ,2)
    elif ch>=61  and ch<=119:return np.arange(61, 119 + 1 ,2)
    elif ch>=121 and ch<=179:return np.arange(121,179 + 1 ,2)
    elif ch>=181 and ch<=239:return np.arange(181,239 + 1 ,2)
    elif ch>=241 and ch<=299:return np.arange(241,299 + 1 ,2)
    elif ch>=301 and ch<=359:return np.arange(301,359 + 1 ,2)
    elif ch>=361 and ch<=419:return np.arange(361,419 + 1 ,2)
    elif ch>=421 and ch<=479:return np.arange(421,479 + 1 ,2)
    else: return -1

# add for MUTES to read two column
def getCollist2(ch):
    if   ch>=1   and ch<=59: return np.arange(1,  119 + 1 ,2)
    elif ch>=61  and ch<=119:return np.arange(1, 119 + 1 ,2)
    elif ch>=121 and ch<=179:return np.arange(121,239 + 1 ,2)
    elif ch>=181 and ch<=239:return np.arange(121,239 + 1 ,2)
    elif ch>=241 and ch<=299:return np.arange(241,359 + 1 ,2)
    elif ch>=301 and ch<=359:return np.arange(241,359 + 1 ,2)
    elif ch>=361 and ch<=419:return np.arange(361,479 + 1 ,2)
    elif ch>=421 and ch<=479:return np.arange(361,479 + 1 ,2)
    else: return -1



def getgrptriglist_inCol(fname):

    f = open(fname, 'w') 

    for one in pos:

        dist = []
        gch = []
        ch = one[0]
        clist = getCollist(ch)        

        for ich in clist:
            d = getDistance(ch,ich)
            dist.append(d)
            gch.append(ich)
            
        dist = np.array(dist)
        gch = np.array(gch)
        
        id = np.argsort(dist)
        dist = dist[id]
        gch = gch[id]

        outstr = str(ch) + "," + str(gch[1]) + "," + str(gch[2]) + ","+ str(gch[3]) + ","+ str(gch[4]) + "\n"
        print " distance of ", ch, " : ", dist[1],",", dist[2], ",", dist[3], ",", dist[4]
        print " nearby ch   ", ch,  " : ", gch[1], ",",  gch[2],   ",", gch[3],  ",", gch[4]
        print outstr
        print ""
        f.write(outstr)


def getgrptriglist_inCol2(fname):

    f = open(fname, 'w') 

    for one in pos:

        dist = []
        gch = []
        ch = one[0]
        clist = getCollist2(ch)        

        for ich in clist:
            d = getDistance(ch,ich)
            dist.append(d)
            gch.append(ich)
            
        dist = np.array(dist)
        gch = np.array(gch)
        
        id = np.argsort(dist)
        dist = dist[id]
        gch = gch[id]

        outstr = str(ch) + "," + str(gch[1]) + "," + str(gch[2]) + ","+ str(gch[3]) + ","+ str(gch[4]) + "\n"
        print " distance of ", ch, " : ", dist[1],",", dist[2], ",", dist[3], ",", dist[4]
        print " nearby ch   ", ch,  " : ", gch[1], ",",  gch[2],   ",", gch[3],  ",", gch[4]
        print outstr
        print ""
        f.write(outstr)



def getColColor(ch):
    if   ch>=1   and ch<=59:  color = cm.hsv(0./8.) 
    elif ch>=61  and ch<=119: color = cm.hsv(1./8.)
    elif ch>=121 and ch<=179: color = cm.hsv(2./8.)
    elif ch>=181 and ch<=239: color = cm.hsv(3./8.)
    elif ch>=241 and ch<=299: color = cm.hsv(4./8.)
    elif ch>=301 and ch<=359: color = cm.hsv(5./8.)
    elif ch>=361 and ch<=419: color = cm.hsv(6./8.)
    elif ch>=421 and ch<=479: color = cm.hsv(7./8.)
    return color

def getDensityColors(resols):
    normal = mpl.colors.Normalize(resols.min(),resols.max())
#    normal = plt.normalize(resols.min(),resols.max())
    colors = plt.cm.jet(normal(resols))
    return colors
        
def getRow(ch):
    p = np.where(pos[:,0]==ch)[0]
    row = pos[p][0][3]
    return row

def getTwrPixAdd(row):
    twr=-1
    pix=-1
    if row>=0 and row<=29:
        twr = twradd[row]
        pix = pixadd[row]
    
    return twr, pix

def getPix(ch):
    #ch->row
    row = getRow(ch)
    pix = pixadd[row]
    return pix

def getTwr(ch):
    #ch->row
    row = getRow(ch)
    twr = twradd[row]
    return twr
    
def getTESPos(ch):
    p = np.where(pos[:,0]==ch)[0]
    x = pos[p][0][1]
    y = pos[p][0][2]
    return x, y


def getDistance(ch1,ch2):
    x1, y1 = getTESPos(ch1)
    x2, y2 = getTESPos(ch2)
    d = sqrt(pow((x1-x2),2) + pow((y1-y2),2))
    return d


def calcTwoPixDistance(col1,pix1,col2,pix2):
    row1 = pixadd.index(pix1)
    row2 = pixadd.index(pix2)
    ch1 = 2*(row1+col1*30)+1
    ch2 = 2*(row2+col2*30)+1
    x1, y1 = getTESPos(ch1)
    x2, y2 = getTESPos(ch2)
    d = sqrt(pow((x1-x2),2) + pow((y1-y2),2))
    return d

def calcDistance(x1,y1,col2,pix2):
    row2 = pixadd.index(pix2)
    ch2 = 2*(row2+col2*30)+1
    x2, y2 = getTESPos(ch2)
    d = sqrt(pow((x1-x2),2) + pow((y1-y2),2))
    return d


def getRectTESPos(ch,facecolor='blue',edgecolor='red',alpha=0.8):
    l = 350.
    x,y = getTESPos(ch)
    rect = plt.Rectangle((x-l/2.,y-l/2.),l,l,fill=True,\
                         facecolor=facecolor,edgecolor=edgecolor,alpha=alpha)
    return rect


def writeTextTESPos(ax,ch,idx,fontsize=8):
    l = 350.
    x,y = getTESPos(ch)
    #    ax.annotate("%d"%idx, xy=(x-l/2.,y-l/2.), xycoords='data',\
    #                xytext=(x-l/2.,y-l/2.), textcoords='data',size=fontsize)
    ax.text(x-l/2.,y-l/2.,"%d"%idx,size=fontsize)


def writeFloatTESPos(ax,ch,idx,fontsize=8, yoffset=0.0, color='k'):
    l = 350.
    x,y = getTESPos(ch)
    #    ax.annotate("%d"%idx, xy=(x-l/2.,y-l/2.), xycoords='data',\
    #                xytext=(x-l/2.,y-l/2.), textcoords='data',size=fontsize)
    ax.text(x-l/2.,y-l/2.+yoffset,"%.1f"%idx,size=fontsize,color=color)

def getPlotRectPlane(i=1,j=1,k=1):
    l = 350.
    #plt.figure()
    ax = plt.subplot(i,j,k, aspect='equal',xlim=(-5000,5000),ylim=(-5000,5000))
    for ch in xrange(1,480,2):
        x,y = getTESPos(ch)
        edgecolor = getColColor(ch)
        rect = plt.Rectangle((x-l/2.,y-l/2.),l,l,fill=False,edgecolor=edgecolor,alpha=0.4)
        ax.add_patch(rect)

    return ax



def getPlotDensityPlane(i=1,j=1,k=1,chans=[1],resols=[1.]):
    l = 350.
    #plt.figure()
    ax = plt.subplot(i,j,k, aspect='equal',xlim=(-5000,5000),ylim=(-5000,5000))
    colors = getDensityColors(resols)
    for i in xrange(0,len(chans),1):
        x,y = getTESPos(chans[i])
        c = colors[i]
        rect = plt.Rectangle((x-l/2.,y-l/2.),l,l,fill=True,color=c,alpha=0.2)
        ax.add_patch(rect)

    return ax


'''
def getRectTESPos(ch):
    l = 350.
    x,y = getTESPos(ch)
    rect = matplotlib.patches.Rectangle((x-l/2.,y-l/2.),l,l)
    return rect

def plotRectPlane():
    l = 350.
    plt.figure()
    ax = plt.subplot(1,1,1, aspect='equal')
    for ch in xrange(1,480,2):
        x,y = getTESPos(ch)
        rect = matplotlib.patches.Rectangle((x-l/2.,y-l/2.),l,l,fill=False)
        ax.add_patch(rect)

    pyplt.xlim(-5000,5000)
    pyplt.ylim(-5000,5000)
    pyplt.grid()
    pyplt.show()
'''
