#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

from __future__ import division
from os import walk
from os import path
from os import mkdir

import numpy

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

tag = "_25M_v2_joined_best/"

CN_input = "CN_joint_log2rs" + tag
BAF_input = "BAF_filtered_data_25M_15/"
validation_output = "segmentation_validation" + tag

if not(path.isdir(validation_output + "/")):
    mkdir(validation_output + "/")

CNlist = []
for (_, _, f) in walk(CN_input):
    CNlist += f

BAFlist = []
for (_, _, f) in walk(BAF_input):
    BAFlist += f

#runall = open(expands_output + "runall.bat", "w")
segments = {}
for CNfile in CNlist:
    (id, sample, _) = CNfile.split("_")
#    if id != "71":
#        continue
    BAFname = id + "_" + sample + "_BAF.txt"
    BAFname = BAFname.translate(None, 'b')
    if not(BAFname in BAFlist):
        print "Couldn't find expected BAF file", BAFname, "from CN file", CNfile
        continue
    if not(id in segments):
        segments[id] = {}
    segments[id][sample] = {}
    cnf = open(CN_input + CNfile, "r")
    for line in cnf:
        if (line.find("chr") != -1):
            continue
        (chr, start, end, rawA, rawB, intA, intB, log2r, nSNPS) = line.split()
        if (chr == "23"):
            continue
        if (chr == "24"):
            continue
        if (end == "inf"):
            end = lsl.getChromosomeMax(int(chr))
        if (intA == intB):
            continue #Double delete; wt; balanced gain: all won't have differential BAFs
        if not(chr in segments[id][sample]):
            segments[id][sample][chr] = []
        segments[id][sample][chr].append([int(start), int(end)])

failfile = open(validation_output + "failures.txt", "w")
summaryfile = open(validation_output + "summary.txt", "w")
failfile.write("patient\tchr\tstart\tend\tsample\tBAFs-.5\n")
summaryfile.write("patient\tsample\tmatches\tantimatches\tfails\n")

for id in segments:
#    if id != "71":
#        continue
    BAFs_by_sample = {}
    BAF_averages = {}
    for sample in segments[id]:
        #if sample != "18992":
        #    continue
        testfile = open("valseg_pjbs.txt", "w")
        if not(sample in BAFs_by_sample):
            BAFs_by_sample[sample] = {}
        if not(sample in BAF_averages):
            BAF_averages[sample] = {}
        for chr in segments[id][sample]:
            for seg in segments[id][sample][chr]:
                segname = (chr, seg[0], seg[1])
                if not(segname in BAFs_by_sample[sample]):
                    BAFs_by_sample[sample][segname] = {}
        BAFname = id + "_" + sample + "_BAF.txt"
        BAFname = BAFname.translate(None, 'b')
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
            pos = int(pos)
            if chr not in segments[id][sample]:
                continue
            for seg in segments[id][sample][chr]:
                #These segments come from ASCAT, which is a both-end-inclusive caller (i.e. calls "5-10", "11-24", etc.)
                if pos < seg[0]:
                    continue
                if pos > seg[1]:
                    continue
                segname = (chr, seg[0], seg[1])
                pos = str(pos)
                if pos in BAFs_by_sample[sample][segname]:
                    pos = pos = "_b"
                BAFs_by_sample[sample][segname][pos] = baf
                if chr=="9" and seg[0] == 23915711:
                    testfile .write(chr + "\t" + pos + "\t" + str(baf) + "\n")
                break
        for segname in BAFs_by_sample[sample]:
            allbafs = []
            for pos in BAFs_by_sample[sample][segname]:
                baf = BAFs_by_sample[sample][segname][pos]
                if (baf<0.5):
                    baf = 1-baf
                allbafs.append(baf)
            if len(allbafs) > 0:
                BAF_averages[sample][segname] = numpy.average(allbafs)
            else:
                BAF_averages[sample][segname] = 0.5 #All wt entries are skipped, like this should be.
    print "Processing patient", id
    lsl.validateSegments(BAFs_by_sample, BAF_averages, validation_output, id, failfile, summaryfile)
        
    
failfile.close()
testfile.close()
summaryfile.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        