#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 12:50:14 2017

@author: lpsmith
"""

from __future__ import division
from os import walk

import lucianSNPLibrary as lsl
import numpy

(seglengths, avgseg) = lsl.getNumSegmentsPerSample()

directories = ["expands_full_input/results/"]
for n in range(100):
    directories.append("expands_full_input_shuffled_" + str(n) + "/results/")

print "average, stdev, median, max, min"
for directory in directories:
    flist = []
    for (_, _, f) in walk(directory):
        flist += f

    lownoisevals = []
    highnoisevals = []
    for f in flist:
        if f.find(".spstats") == -1:
            continue
        (patient, sample, tag) = f.split("_")
        statfile = open(directory + f, "r")
        for line in statfile:
            if line.find("expands") != -1:
                continue
            if line.find("Mean") != -1:
                continue
            splitline = line.split()
            if len(splitline) < 2:
                continue
            if seglengths[(patient, sample)] < avgseg:
                lownoisevals.append(float(splitline[1]))
            else:
                highnoisevals.append(float(splitline[1]))
        statfile.close()
    if (len(highnoisevals) > 0):
        print directory, "High:", numpy.average(highnoisevals), numpy.std(highnoisevals), numpy.median(highnoisevals), numpy.max(highnoisevals), numpy.min(highnoisevals)
        print directory, "Low:", numpy.average(lownoisevals), numpy.std(lownoisevals), numpy.median(lownoisevals), numpy.max(lownoisevals), numpy.min(lownoisevals)
        
        lsl.createPrintAndSaveHistogram(highnoisevals, "highnoiseout.txt", 0.01, xdata="noise")
        lsl.createPrintAndSaveHistogram(lownoisevals, "lownoiseout.txt", 0.01, xdata="noise")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        