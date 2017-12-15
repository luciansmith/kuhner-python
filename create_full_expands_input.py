#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

from __future__ import division
from os import walk

import math
import numpy
from scipy.interpolate import interp1d
import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not
rejoined = ""
#rejoined = "_rejoined"

points = [0.1, 1.0,2.0,3.0,4.0,5.0]
logR = [-2.5, -0.44, 0.0, 0.163, 0.282, 0.37]
# -0.44 from peak of CN_rejoined_histograms/loss_21-100000000.txt
# Other values from separate_histograms.py running on CN_rejoined_histograms/balanced_gain_hist_21-100000000.txt
values = []
maxval = 0
for v in logR:
  values.append(2.0 * math.pow(2.0,v))
f2 = interp1d(points,values,kind="cubic")

CN_input = "CN_calc_log2rs" + rejoined + "/"
BAF_input = "BAF_filtered_data_15/"
expands_output = "expands_full_input" + rejoined + "/"

CNlist = []
for (_, _, f) in walk(CN_input):
    CNlist += f

BAFlist = []
for (_, _, f) in walk(BAF_input):
    BAFlist += f

#runall = open(expands_output + "runall.bat", "w")
for CNfile in CNlist:
    (id, sample, _) = CNfile.split("_")
    if id != "99":
        continue
    BAFname = id + "_" + sample + "_BAF.txt"
    if not(BAFname in BAFlist):
        print "Couldn't find expected BAF file", BAFname, "from CN file", CNfile
        continue
    cnf = open(CN_input + CNfile, "r")
    print "Writing", id+"_"+sample+"_CN.txt"
    CN_infile = open(expands_output + id + "_" + sample + "_CN.txt", "w")
    CN_infile.write("chr\tstartpos\tendpos\tCN_Estimate\n")
    segments = {}
    for line in cnf:
        if (line.find("chr") != -1):
            continue
        (chr, start, end, Xl2r, call, nl2r, log2r,stdev) = line.split()
        if (chr == "23"):
            continue
        if (chr == "24"):
            continue
        if (end == "inf"):
            end = lsl.getChromosomeMax(int(chr))
        if not(chr in segments):
            segments[chr] = []
        segments[chr].append([int(start), int(end), list()])
        try:
            log2r = float(log2r)
        except:
            continue
        if (log2r > logR[len(logR)-1]):
            continue
        val = 2.0 * math.pow(2.0,log2r)
        # interpolation
        val = f2(val)
        CN_infile.write(chr + "\t" + start + "\t" + str(end) + "\t" + str(val) + "\n")

    CN_infile.close()
    baffile = open(BAF_input + BAFname, "r")
    for line in baffile:
        if (line.find("BAF") != -1):
            continue
        (snpid, chr, pos, baf) = line.split()
        if (baf=="?"):
            continue
        if (chr == "23"):
            continue
        if (chr == "24"):
            continue
        baf = float(baf)
        if (baf < 0.5):
            baf = 1-baf
        pos = int(pos)
        for segment in segments[chr]:
            if pos <= segment[0]:
                continue
            if pos > segment[1]:
                continue
            segment[2].append(baf)
            break
    print "Writing", id+"_"+sample+"_BAF.txt"
    BAF_infile = open(expands_output + id + "_" + sample + "_BAF.txt", "w")
    BAF_infile.write("chr\tstartpos\tendpos\tAF_Tumor\tPN_B\n")
    for c in range(1,22):
        for segment in segments[str(c)]:
            if len(segment[2]) == 0:
                continue
            BAF_infile.write(str(c) + "\t")
            BAF_infile.write(str(segment[0]) + "\t")
            BAF_infile.write(str(segment[1]) + "\t")
            BAF_infile.write(str(numpy.average(segment[2])) + "\t1\n")
#            if (c == 3 and segment[0] == 60384729):
#                print "BAF values for chr3,", segment[0], ",", segment[1]
#                print segment[2]
    BAF_infile.close()
    Rout   = open(expands_output + id + "_" + sample + "_analyze.R", "w")
    Rout.write("library('expands')\n")
    Rout.write("runExPANdS('" + id + "_" + sample + "_BAF.txt','" + id + "_" + sample + "_CN.txt', snvF='results" + rejoined + "/" + id + "_" + sample + "_SPs.txt')")
    Rout.close()
    Rout   = open(expands_output + "run_" + id + "_" + sample + "_analyze.sge", "w")
    Rout.write("module load modules modules-init modules-gs gmp/5.0.2 mpfr/latest mpc/0.8.2 gcc/latest R/latest java_jdk/latest\n")
    Rout.write("cd R/" + expands_output + "\n")
    Rout.write("R CMD BATCH " + id + "_" + sample + "_analyze.R\n")
    Rout.close()
    #runall.write("qsub run_" + id + "_" + sample + "_analyze.sge\n")

#runall.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        