#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:02:52 2016

@author: lpsmith
"""

#Create parallel BAF segmentation files to match the CN segmentation files

from __future__ import division
from os import walk
import sys

import lucianSNPLibrary as lsl
import numpy

directories = ["expands_full_input_alldupes_noLOH/results/"]#, "shuffled0/results/"]
for n in range(100):
    directories.append("shuffled_output/expands_full_input_shuffled_noLOH_" + str(n) + "/results/")

outfile = open("spstats.summary.txt", "w")
outfile.write("Run\taverage\tstdev\tmedian\tmax\tmin\n")

for directory in directories:
    flist = []
    for (_, _, f) in walk(directory):
        flist += f
    
    print "Processing", directory
    noisevals = []
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
            noisevals.append(float(splitline[1]))
        statfile.close()
    if (len(noisevals) > 0):
        outfile.write(directory + "\t")
        outfile.write(str(numpy.average(noisevals)) + "\t")
        outfile.write(str(numpy.std(noisevals)) + "\t")
        outfile.write(str(numpy.median(noisevals)) + "\t")
        outfile.write(str(numpy.max(noisevals)) + "\t")
        outfile.write(str(numpy.min(noisevals)) + "\n")
        lsl.createPrintAndSaveHistogram(noisevals, "noiseout.txt", 0.01, xdata="noise")
outfile.close()
of = open("noisevals.txt", "w")
of.write(str(noisevals))
of.close()



















