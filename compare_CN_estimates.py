# convert Patty's data to Expands input
# using Xiaohong interpolation 
# this version also writes the Expands R script
# one sample per run (unlike Pierre data)

import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
import numpy as np
import sys

# produce interpolation curve (info provided by Xiaohong)
points = [1.0,2.0,3.0,4.0]
logR = [-0.55, 0.0, 0.395, 0.53]
values = []
for v in logR:
  values.append(2.0 * math.pow(2.0,v))
f2 = interp1d(points,values,kind="cubic")

# read a Patty logr file
infilename = "1047_18966_partekSegment.txt"
partekfile = open(infilename,"r")

# set up the R script
comparefile = open("CN_compare.txt","w")
comparefile.write("chrom\tstartpos\tendpos\tCN_partek\tCN_kenji\n")

#read the kenji estimates
kenjifile = open("combined_results7.merged_calls", "r")
kenjidata = kenjifile.readlines()
cns_for_1047 = []
cns_for_1047
for line in kenjidata:
    line = line.rstrip().split()
    if (len(line) == 6):
        if (line[0][0:4] == "1047"):
            cns_for_1047.append([line[1], int(line[2]), int(line[3]), int(line[4]), int(line[5])])

partekdata = partekfile.readlines()
partekfile.close()

for line in partekdata:
    line = line.rstrip().split()
    (chr,startpos,endpos,id,classification,length,mean,nummarkers,pval) = line[0:9]
    chr = chr[3:len(chr)]
    if (chr=="X"):
        chr = "23"
    if (chr=="Y"):
        chr = "24"
    outline = chr + "\t" + startpos + "\t" + endpos
    startpos = int(startpos)
    endpos = int(endpos)
    val = float(mean)
    val = 2.0 * math.pow(2.0,val)
# interpolation
    if (1.0 < val < 4.0):
        cn = f2(val)
        outline += "\t" + str(cn)
        #print "logr",mean,"estimate",cn,classification
    else:
        print "logr",mean,"failed",classification
        continue
    #Now find that value in the kenji data
    kfinal = "??"
    print "Looking for overlap, chromosome ", chr, ", ",startpos, "-", endpos
    for kline in cns_for_1047:
        if (kline[0] != chr):
            continue
        (kstart, kend, ksize, kcall) = kline[1:5]
        if ((startpos < kstart or startpos > kend) and \
            (endpos < kstart or endpos > kend)):
            if (startpos < kstart and endpos > kend):
                print "complete overlap!"
            else:
                print "no overlap: ",kstart,"-",kend
                continue
        print "Overlap found: ", kstart, "-", kend
        if (endpos > kend):
            print "unequal break\n"
        if (kcall == 0):
            kfinal = "0"
        elif(kcall == 1):
            kfinal = "1"
        elif(kcall == 2):
            kfinal = "2"
        elif(kcall == 3):
            kfinal = "2"
        elif(kcall == 4):
            kfinal = "3"
        elif(kcall == 5):
            kfinal = "4"
        elif(kcall == 6):
            kfinal = "6"
        else:
            kfinal = "unknown"
        break
    outline += "\t" + kfinal + "\n"
    comparefile.write(outline)
comparefile.close()
