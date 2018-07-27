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


def validateSegments(BAFs_by_sample, validation_output, id, failfile, summaryfile):
    #assumes that any double deletion is already removed.
    output = {}
    samples = list(BAFs_by_sample.keys())
    sampairs = set()
    summatches = 0
    sumantimatches = 0
    sumfails = 0
    samplefails = {}
    for sample in samples:
        samplefails[sample] = [0, 0, 0]
    for n in range(0,len(samples)-1):
        sample1 = samples[n]
        for m in range(n+1, len(samples)):
            sample2 = samples[m]
            matches = []
            for segname in BAFs_by_sample[sample1]:
                if not(segname in BAFs_by_sample[sample2]):
                    continue
                match = 0
                antimatch = 0
                for pos in BAFs_by_sample[sample1][segname]:
                    if not(pos in BAFs_by_sample[sample2][segname]):
                        continue
                    baf1 = BAFs_by_sample[sample1][segname][pos]
                    baf2 = BAFs_by_sample[sample2][segname][pos]
#                    if 0.4 < baf1 and baf1 < 0.6:
#                        continue
#                    if 0.4 < baf2 and baf2 < 0.6:
#                        continue
                    if baf1 < 0.5 and baf2 < 0.5:
                        match += 1
                    elif baf1 > 0.5 and baf2 > 0.5:
                        match += 1
                    elif baf1 < 0.5 and baf2 > 0.5:
                        antimatch += 1
                    elif baf1 > 0.5 and baf2 < 0.5:
                        antimatch += 1
                    #print(baf1, baf2, match, antimatch)
                matches.append((match, antimatch))
                nmax = match + antimatch
                if nmax >= 2:
                    if match/nmax > 0.9:
                        samplefails[sample1][0] += 1
                        samplefails[sample2][0] += 1
                    elif antimatch/nmax > 0.9:
                        samplefails[sample1][1] += 1
                        samplefails[sample2][1] += 1
                    else:
                        samplefails[sample1][2] += 1
                        samplefails[sample2][2] += 1
                if not(segname in output):
                    output[segname] = {}
                output[segname][sample1+"_"+sample2] = (match, antimatch)
                sampairs.add(sample1 + "_" + sample2)
    for sample in samples:
        summaryfile.write(id + "\t" + sample + "\t" + str(samplefails[sample][0]) + "\t" + str(samplefails[sample][1]) + "\t" +  str(samplefails[sample][2])  + "\n")
    outfile = open(validation_output + id + "_segvalidate.txt", "w")
    outfile.write("patient\tchr\tstart\tend\tmax_nBAFs\tmatches\tantimatches\tfails")
    sampairs = sorted(sampairs)
    for sampair in sampairs:
        outfile.write("\t" + sampair + " match\t" + sampair + " anti-match")
    outfile.write("\n")
    for segname in output:
        outfile.write(id + "\t" + segname[0] +"\t" + str(segname[1]) +"\t" + str(segname[2]))
        outline = ""
        nummatches = 0
        numantimatches = 0
        numfails = 0
        maxBAFs = 0
        for sampair in sampairs:
            if sampair in output[segname]:
                match = output[segname][sampair][0]
                antimatch = output[segname][sampair][1]
                outline += "\t" + str(match) + "\t" + str(antimatch)
                numBAFs = match+antimatch
                if maxBAFs<numBAFs:
                    maxBAFs = numBAFs
                if numBAFs<=1:
                    continue
                elif (match/numBAFs) > .9:
                    #print("Match", match, antimatch, numBAFs, match/numBAFs)
                    nummatches += 1
                elif (antimatch/numBAFs) > .9:
                    numantimatches += 1
                    #print("Antimatch", match, antimatch, numBAFs, match/numBAFs)
                else:
                    numfails += 1
                    (sample1, sample2) = sampair.split("_")
                    failfile.write(id + "\t" + segname[0] + "\t" + str(segname[1]) +"\t" + str(segname[2]) + "\t" + sample1)
                    for pos in BAFs_by_sample[sample1][segname]:
                        failfile.write("\t" + str(BAFs_by_sample[sample1][segname][pos]))
                    failfile.write("\n")
                    failfile.write(id + "\t" + segname[0] + "\t" + str(segname[1]) +"\t" + str(segname[2]) + "\t" + sample2)
                    for pos in BAFs_by_sample[sample2][segname]:
                        failfile.write("\t" + str(BAFs_by_sample[sample2][segname][pos]))
                    failfile.write("\n")
                    #print("Fail", match, antimatch, numBAFs, match/numBAFs)
            else:
                outline += "\t--\t--"
        outfile.write("\t" + str(maxBAFs) + "\t" + str(nummatches) + "\t" + str(numantimatches) + "\t" + str(numfails) + outline + "\n")
        summatches += nummatches
        sumantimatches += numantimatches
        sumfails += numfails
    outfile.close()
    summaryfile.write(id + "\tall\t" + str(summatches) + "\t" + str(sumantimatches) + "\t" + str(sumfails) + "\n")



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
    validateSegments(BAFs_by_sample, validation_output, patient, failfile, summaryfile)


failfile.close()
summaryfile.close()

















