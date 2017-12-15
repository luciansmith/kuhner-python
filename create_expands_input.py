# -*- coding: utf-8 -*-
"""
Created on Wed May 11 17:58:38 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
import math
from scipy.interpolate import interp1d
import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not
rejoined = ""
#rejoined = "_rejoined"

inputdir = "BAF_persegment" + rejoined + "/"
outdir = "expands" + rejoined + "_input_alldupes/"

points = [0.1, 1.0,2.0,3.0,4.0,5.0]
logR = [-2.5, -0.44, 0.0, 0.163, 0.282, 0.37]
# -0.44 from peak of CN_rejoined_histograms/loss_21-100000000.txt
# Other values from separate_histograms.py running on CN_rejoined_histograms/balanced_gain_hist_21-100000000.txt
values = []
for v in logR:
  values.append(2.0 * math.pow(2.0,v))
f2 = interp1d(points,values,kind="cubic")

flist = []
#SNPfiles.append(["1034", "20008"])
for (_, _, f) in walk(inputdir):
    flist += f

for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 2):
        continue
    patient = split[0]
    sample = split[1].replace(".txt","")
    infile = open(inputdir + f, "r")
    bafout = open(outdir + patient + "_" + sample + "_BAF.txt", "w")
    cnout  = open(outdir + patient + "_" + sample + "_CN.txt",  "w")
    bafout.write("chr\tstartpos\tendpos\tAF_Tumor\tPN_B\n")
    cnout.write ("chr\tstartpos\tendpos\tCN_Estimate\n")
    for line in infile:
        if (line.find("chr") != -1):
            continue
        (chr, start, end, log2r, baf) = line.strip().split()
        if (end=="inf"):
            end = str(lsl.getChromosomeMax(int(chr)))
        if (log2r == "?"):
            continue
        log2r = float(log2r)
        val = 2.0 * math.pow(2.0,log2r)
        # interpolation
        val = f2(val)
        if (baf != "---"):
            bafout.write(chr + "\t" + start + "\t" + end + "\t" + baf + "\t" + "1\n")
            cnout.write(chr + "\t" + start + "\t" + end + "\t" + str(val) + "\n")
    bafout.close()
    cnout.close()
    Rout   = open(outdir + patient + "_" + sample + "_analyze.R", "w")
    Rout.write("library('expands')\n")
    Rout.write("runExPANdS('" + patient + "_" + sample + "_BAF.txt','" + patient + "_" + sample + "_CN.txt', snvF='results" + rejoined + "/" + patient + "_" + sample + "_SPs.txt')")
    Rout.close()
    Rout   = open(outdir + "run_" + patient + "_" + sample + "_analyze.sge", "w")
    Rout.write("module load modules modules-init modules-gs gmp/5.0.2 mpfr/latest mpc/0.8.2 gcc/latest R/latest java_jdk/latest\n")
    Rout.write("cd R/" + outdir + "\n")
    Rout.write("R CMD BATCH " + patient + "_" + sample + "_analyze.R\n")
    Rout.close()
    
    
    
    