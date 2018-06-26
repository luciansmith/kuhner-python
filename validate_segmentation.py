#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Cproreated on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

from __future__ import division
from os import walk
from os import path
from os import mkdir
from os.path import isfile

import numpy

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

tag = "_g500_better_ploidy/"

BAF_input = ["BAF_filtered_data_25M_only_40_65/", "BAF_filtered_data_1M_only_40_65/", "BAF_filtered_data_Pilot_40_65/"]
validation_output = "segmentation_validation" + tag

use_nonints = True

if use_nonints:
    CN_input = "noninteger_processed_CNs/"
else:
    CN_input = "CN_joint_log2rs" + tag

if not(path.isdir(validation_output + "/")):
    mkdir(validation_output + "/")

CNlist = []
for (_, _, f) in walk(CN_input):
    CNlist += f

BAFlist = []
for bafin in BAF_input:
    for (_, _, f) in walk(bafin):
        BAFlist += f

#runall = open(expands_output + "runall.bat", "w")
segments = {}
for CNfile in CNlist:
    if use_nonints:
        (patient, sample, gamma, ploidy, __, __) = CNfile.split("_")
        if ploidy != lsl.getBestPloidyFor(patient, sample):
            continue
    else:
        (patient, sample, __) = CNfile.split("_")
#    if patient != "71":
#        continue
    BAFname = patient + "_" + sample + "_BAF.txt"
    BAFname = BAFname.replace('b', '')
    if not(BAFname in BAFlist):
        print("Couldn't find expected BAF file", BAFname, "from CN file", CNfile)
        continue
    if not(patient in segments):
        segments[patient] = {}
    segments[patient][sample] = {}
    cnf = open(CN_input + CNfile, "r")
    for line in cnf:
        if (line.find("chr") != -1):
            continue
        if use_nonints:
            (__, __, chr, start, end, __, __, intA, intB) = line.rstrip().split()
        else:
            (chr, start, end, __, __, intA, intB, __, __) = line.rstrip().split()
        if (chr == "23"):
            continue
        if (chr == "24"):
            continue
        if (end == "inf"):
            end = lsl.getChromosomeMax(int(chr))
        if (intA == intB):
            continue #Double delete; wt; balanced gain: all won't have differential BAFs
        if not(chr in segments[patient][sample]):
            segments[patient][sample][chr] = []
        segments[patient][sample][chr].append([int(start), int(end)])

failfile = open(validation_output + "failures.txt", "w")
summaryfile = open(validation_output + "summary.txt", "w")
failfile.write("patient\tchr\tstart\tend\tsample\tBAFs-.5\n")
summaryfile.write("patient\tsample\tmatches\tantimatches\tfails\n")

for patient in segments:
#    if patient != "71":
#        continue
    BAFs_by_sample = {}
    for sample in segments[patient]:
        #if sample != "18992":
        #    continue
        if not(sample in BAFs_by_sample):
            BAFs_by_sample[sample] = {}
        for chr in segments[patient][sample]:
            for seg in segments[patient][sample][chr]:
                segname = (chr, seg[0], seg[1])
                if not(segname in BAFs_by_sample[sample]):
                    BAFs_by_sample[sample][segname] = {}
        BAFname = patient + "_" + sample + "_BAF.txt"
        BAFname = BAFname.replace('b', '')
        baffilename = BAF_input[0] + BAFname
        if not isfile(baffilename):
            baffilename = BAF_input[1] + BAFname
        if not isfile(baffilename):
            baffilename = BAF_input[2] + BAFname
        baffile = open(baffilename, "r")
        for line in baffile:
            if (line.find("BAF") != -1):
                continue
            (snpid, chr, pos, baf) = line.rstrip().split()
            if (baf=="?"):
                continue
            if (chr == "23"):
                continue
            if (chr == "24"):
                continue
            baf = float(baf)
            pos = int(pos)
            if chr not in segments[patient][sample]:
                continue
            for seg in segments[patient][sample][chr]:
                #These segments come from ASCAT, which is a both-end-inclusive caller (i.e. calls "5-10", "11-24", etc.)
                if pos < seg[0]:
                    continue
                if pos > seg[1]:
                    continue
                segname = (chr, seg[0], seg[1])
                pos = str(pos)
                if pos in BAFs_by_sample[sample][segname]:
                    pos = pos + "_b"
                BAFs_by_sample[sample][segname][pos] = baf
                break
    print("Processing patient", patient)
    lsl.validateSegments(BAFs_by_sample, validation_output, patient, failfile, summaryfile)


failfile.close()
summaryfile.close()

















