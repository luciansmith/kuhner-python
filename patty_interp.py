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

pid = sys.argv[1]

# read a Patty logr file
infilename = pid + "_partekSegment.txt"
infile = open(infilename,"r")

# set up the R script
rfile = open("expandsruns.R","w")
rfile.write("library('expands')\n")

indata = infile.readlines()
infile.close()

outfilename = "data/" + pid + "_CN.txt"
outfile = open(outfilename,"w")
outline = "chr\tstartpos\tendpos\tCN_Estimate\n"
outfile.write(outline)
numfail = 0
for line in indata:
  line = line.rstrip().split()
  (chr,startpos,endpos,id,classification,length,mean,nummarkers,pval) = line[0:9]
  outline = chr + "\t" + startpos + "\t" + endpos + "\t" 
  val = float(mean)
  val = 2.0 * math.pow(2.0,val)
# interpolation
  if (1.0 < val < 4.0):
    cn = f2(val)
    outline += str(cn) + "\n"
    outfile.write(outline)
    print "logr",mean,"estimate",cn,classification
  else:
    numfail += 1
    print "logr",mean,"failed",classification
    continue
outfile.close()

# report on success
print "Processed",len(indata),"segments with",numfail,"failures"

# write the R script line to run this file and its partner
bafname = "data/" + pid + "_BAF.txt"
spname = "results/"+ pid + "_SPs.txt"
rline = "runExPANdS('" + bafname + "','" + outfilename
rline += "', snvF='" + spname + "')\n"
rfile.write(rline)
