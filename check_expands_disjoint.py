#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 12:29:30 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
import sys

import lucianSNPLibrary as lsl
import numpy

directories = []#["expands_full_input_alldupes_noLOH/results/", "expands_full_input/results/"] #, "shuffled0/results/"]
for n in range(100):
    directories.append("shuffled_output/expands_full_input_shuffled_noLOH_" + str(n) + "/results/")

outfile = open("sps.disjoint.summary.txt", "w")
outfile.write("Directory\tPatient\tSample\tNumDisjoints\tNumPossibleDisjoints\n")

def addDisjoints(low, high, disjoints, asnts):
    for asnt in asnts:
        if asnt[high] == "1":
            disjoints[1] += 1
            if asnt[low] != "1":
                disjoints[0] += 1

overall = [0, 0, 0, 0, 0]
for directory in directories:
    flist = []
    for (_, _, f) in walk(directory):
        flist += f
    
    print "Processing", directory
    maxpvals = []
    logmaxpvals = []
    for f in flist:
        if f.find(".spstats") != -1:
            continue
        if f.find(".sps.cbs") != -1:
            continue
        if f.find(".sps") == -1:
            continue
        (patient, sample, __) = f.split("_")
        spsfile = open(directory + f, "r")
        spsfile.readline() #'Expands version 1.7.2' (or whatever)
        headers = spsfile.readline().rstrip().split()
        if len(headers) < 16:
            continue
        first = 15
        last = len(headers)-1
        SPs = headers[first:last]
        for n in range(len(SPs)):
            SPs[n] = float(SPs[n][3:len(SPs[n])])
        asnts = []
        for line in spsfile:
            splitline = line.rstrip().split()
            asnts.append(splitline[first:last])
        spsfile.close()
        disjoints = [0, 0]
        for n in range(len(SPs)):
            for m in range(n+1, len(SPs)):
                sp1 = SPs[n]
                sp2 = SPs[m]
                if (sp1+sp2 <= 1):
                    continue
                if sp1 < sp2:
                    addDisjoints(n, m, disjoints, asnts)
                else:
                    addDisjoints(m, n, disjoints, asnts)
        outfile.write(directory + "\t" + patient + "\t" + sample + "\t" + str(disjoints[0]) + "\t" + str(disjoints[1]) + "\n")
        if disjoints[1] == 0:
            overall[0] += 1 #untestable
        elif disjoints[0] == disjoints[1]:
            overall[1] += 1 #always failed
        elif disjoints[0] == 0:
            overall[4] += 1 #always passed
        elif disjoints[0]/disjoints[1] > .5:
            overall[2] += 1 #failed more than passed
        else:
            overall[3] += 1 #passed more than failed
outfile.close()
print overall
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        