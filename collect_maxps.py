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

directories = ["expands_full_input_alldupes_noLOH/results/"] #, "shuffled0/results/"]
for n in range(100):
    directories.append("shuffled_output/expands_full_input_shuffled_noLOH_" + str(n) + "/results/")

outfile = open("sps.maxP.summary.txt", "w")
outfile.write("Run\taverage\tstdev\tmedian\tmax\tmin\taverage-log\tstdev-log\tmedian-log\tmax-log\tmin-log\n")

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
        for line in spsfile:
            if line.find("expands") != -1:
                continue
            if line.find("chr") != -1:
                continue
            splitline = line.split()
            if len(splitline) < 9:
                continue
            maxp = splitline[8]
            try:
                maxp = float(maxp)
            except:
                continue
            maxpvals.append(maxp)
            logmaxpvals.append(numpy.log(maxp))
        spsfile.close()
    if (len(maxpvals) > 0):
        outfile.write(directory + "\t")
        outfile.write(str(numpy.average(maxpvals)) + "\t")
        outfile.write(str(numpy.std(maxpvals)) + "\t")
        outfile.write(str(numpy.median(maxpvals)) + "\t")
        outfile.write(str(numpy.max(maxpvals)) + "\t")
        outfile.write(str(numpy.min(maxpvals)) + "\t")
        outfile.write(str(numpy.average(logmaxpvals)) + "\t")
        outfile.write(str(numpy.std(logmaxpvals)) + "\t")
        outfile.write(str(numpy.median(logmaxpvals)) + "\t")
        outfile.write(str(numpy.max(logmaxpvals)) + "\t")
        outfile.write(str(numpy.min(logmaxpvals)) + "\n")
        #lsl.createPrintAndSaveHistogram(maxpvals, "maxpvals.txt", 1, xdata="MaxP")
        #lsl.createPrintAndSaveHistogram(logmaxpvals, "logmaxpvals.txt", .1, xdata="logMaxP")
outfile.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        