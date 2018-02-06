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
from os import readlink
from os.path import isfile
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

BAF_links = False
if BAF_links:
    BAF_dir = "gamma_template/"
else:
    BAF_dir = "pASCAT_input_combined_all/"

subdirs = ["diploid", "tetraploid"]
#subdirs = ["diploid"]
subdirdict = {}
for subdir in subdirs:
    subdirdict[subdir] = []
gamma_outdir = "gamma_test_output/"
#gamma_outdir = "gamma_test_output/"
outdir = "gamma_test_output/analysis_compare/"
#gamma_list = ["100", "150", "200", "250", "300", "350", "400", "450", "500", "600"]
#gamma_list = ["50"]
#gamma_list = ["Test"]
#gamma_list = ["100", "500", "1000", "3000"]
gamma_list = ["0", "50", "100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000"]
#gamma_list = []
bafrawdata = {}
patient_samples = {}

onlyonepatient = False
onepatient = "163"

#twopatients = [("521", "252")]
#samples from patient 360:
onepatientsamples = ["360_21620", "360_18370", "360_18354", "360_18358", "360_21608", "360_21606", "360_21604", "360_21602", "360_18356", "360_21600", "360_24412", "360_18350", "360_18444", "360_18446", "360_18352", "360_18366", "360_21618", "360_21598", "360_18362", "360_21596", "360_21614", "360_21594", "360_21616", "360_21592", "360_21610", "360_21612", "360_24409", "360_24810", "360_24813", "360_18372"]
ignore_cnvis = True
mirror_percentages = False

#median for patient 521:  .5154885
#Reasonable cutoffs for 521: 0.48 and 0.563 (about 750k left each)
#A more agressive but balanced cutoff:  0.45 and 0.6  (about 215k left each)
#The best cutoff is no cutoff at all:  0.5 and 0.5
bafErrorBarLow = 0.5
bafErrorBarHigh = 0.5

print("Using BAF filtering for comparison from ", str(bafErrorBarLow), "to", str(bafErrorBarHigh))

bafWtLow = 0.4
bafWtHigh = 0.65

#original cutoff: 0.35 to 0.65, or 0.15 out from 0.5.  The above values account for dye bias.

def readBafNormal(patient):
    bafrawdata = {}
    bafwt = {}
#    allcnvis = []
    if BAF_links:
        bafnormal = readlink(BAF_dir + patient + "_Normal_BAF.txt")
    else:
        bafnormal = BAF_dir + patient + "_Normal_BAF.txt"
    if not(isfile(bafnormal)):
        print("ERROR:  no Normal BAF file found for patient", patient)
        return ({}, {})
    bafnormal = open(bafnormal, "r")
    print("Reading BAF normal data for patient", patient)
    for line in bafnormal:
        lvec = line.split()
        if line.find("Chr") != -1:
            continue
        if ignore_cnvis and line.find("cnvi") != -1:
            continue
        try:
            value = float(lvec[3])
        except:
            continue
#        if (line.find("cnvi") != -1):
#            allcnvis.append(value)
        if (value < bafWtLow or value > bafWtHigh):
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if chr not in bafrawdata:
            bafrawdata[chr] = {}
        if chr not in bafwt:
            bafwt[chr] = {}
        bafrawdata[chr][pos] = {}
        bafwt[chr][pos] = value
    bafnormal.close()
#    print("All normal CNVI bafs:")
#    lsl.createPrintAndSaveHistogram(allcnvis, "", .01)
    return bafrawdata, bafwt

def readBafSamples(baffile, bafrawdata):
    labels = []
    all_samples = []
    if BAF_links:
        baffile = readlink(BAF_dir + patient + "_BAF.txt")
    else:
        baffile = BAF_dir + patient + "_BAF.txt"
    if not(isfile(baffile)):
        print("ERROR:  no BAF file found for patient", patient)
        bafrawdata = {}
        return
    print("Reading BAF sample data for patient", patient)
    baffile = open(baffile, "r")
    allbafs = {}
    allbafs["cnvi"] = []
    allbafs["normal"] = []
    for line in baffile:
        lvec = line.split()
        if line.find("Chr") != -1:
            labels = lvec
            for l in range(2,len(labels)):
                all_samples.append(labels[l].split('"')[1])
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if chr not in bafrawdata:
            continue
        if pos not in bafrawdata[chr]:
            continue
        for p in range(3,len(lvec)):
            sample = labels[p-1].split('"')[1]
            try:
                bafrawdata[chr][pos][sample] = float(lvec[p])
                if line.find("cnvi") != -1:
                    allbafs["cnvi"].append(float(lvec[p]))
                else:
                    bafrawdata[chr][pos][sample] = float(lvec[p])
                    allbafs["normal"].append(float(lvec[p]))
            except:
                continue
    baffile.close()
    return allbafs, all_samples

def getIsegsFromCopynumberFileFor(patient):
    #This function reads the lowest-gamma-value copynumber file it can, reads it, and returns it as 'isegs'.
    isegs = {}
    glist = ("100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000")
    gindex = 0
    gamma = glist[gindex]

    root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"
    isegfilename = root_dir + patient + "_copynumber_segments.txt"
    while not(isfile(isegfilename)) and gindex<len(glist):
        gamma = glist[gindex]
        root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"
        isegfilename = root_dir + patient + "_copynumber_segments.txt"
        gindex += 1
    if gindex >= len(glist):
        print("Cannot find any copynumber file for patient", patient)
        return {}
    isegfile = open(isegfilename, "r")
    for line in isegfile:
        if line.find("Chr") != -1:
            continue
        (chr, start, end, nlogr, nbaf) = line.split()
        nbaf = int(nbaf)
        if nbaf <10:
            continue
        start = int(start)
        end = int(end)
        if chr not in isegs:
            isegs[chr] = []
        isegs[chr].append([start, end])
    return isegs

def readSegmentationFile(ascsegfile, all_samples):
    totsca = {}
    totsca["overall"] = set()
    for sample in all_samples:
        totsca[sample] = set()
    osegs = {}
    for line in open(ascsegfile):
        if line.find("Chr") != -1:
            continue
        (__, sample, chr, start, end, nprobes, mbaf, logr, nA, nB) = line.split()
        sample = sample.split('"')[1]
        start = int(start)
        end = int(end)
        if (nA == nB):
            continue
        if sample not in osegs:
            osegs[sample] = {}
        if chr not in osegs[sample]:
            osegs[sample][chr] = []
        osegs[sample][chr].append((start, end))
        totsca[sample].add((chr, start, end))
        totsca["overall"].add((chr, start, end))
    return (osegs, totsca)

def addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, call):
    if chr == "":
        return
    if chr == "23":
        return
    if chr == "24":
        return
    if chr=="0":
        #print("There's a chromosome zero segment, weirdly:", full_sample, chr, start, end)
        return
    start = int(start)
    end = int(end)
    added = False
    svec = full_sample.split("-")
    if len(svec) < 4:
        svec = full_sample.split("_")
    patient = svec[0]
    sample = svec[1]
    if patient not in Xiaohong_segments:
        Xiaohong_segments[patient] = {}
    if sample not in Xiaohong_segments[patient]:
        Xiaohong_segments[patient][sample] = {}
    if chr not in Xiaohong_segments[patient][sample]:
        Xiaohong_segments[patient][sample][chr] = []
    else:
        lastseg = Xiaohong_segments[patient][sample][chr][-1]
        if lastseg[2] == call and start - lastseg[1] < 5000:
            lastseg[1] = end
            added = True
    if not added:
        Xiaohong_segments[patient][sample][chr].append([start, end, call])

def readXiaohongWGSLOHFile(f, Xiaohong_segments):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (__, full_sample, __, chr, start, end) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, "CNLOH")
    xfile.close()

def readXiaohong1MLOHFile(f, Xiaohong_segments):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (full_sample, chr, start, end, __, __) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, "CNLOH")
    xfile.close()

def readXiaohongCopynumFile(f, Xiaohong_segments):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        lvec = line.split()
        if len(lvec) == 7:
            (full_sample, chr, start, end, __, __, call) = line.split()
        elif len(lvec) == 10:
            (__, full_sample, __, chr, start, end, __, __, __, call) = line.split()
        else:
            print("Incorrect line length:", line)
            continue
        if call=="41":
            call = "Double_d"
        elif call=="34":
            call = "Balanced_gain"
        elif call=="22" or call=="23":
            call = "Loss"
        elif call=="32" or call=="33":
            call = "Gain"
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, call)
    xfile.close()

def readAllXiaohongSegmentation():
    print("Reading all Xiaohong segmentation.")
    Xiaohong_segments = {}
    files = []
    Xdir_WGS = "Xiaohong_WGS_segmentation/"
    Xdir_1M = "CN_Xiaohong_segmentation/"
    for (__, __, f) in walk(Xdir_WGS):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("test") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohongWGSLOHFile(Xdir_WGS + f, Xiaohong_segments)
        else:
            readXiaohongCopynumFile(Xdir_WGS + f, Xiaohong_segments)
    files = []
    for (__, __, f) in walk(Xdir_1M):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohong1MLOHFile(Xdir_1M + f, Xiaohong_segments)
        else:
            readXiaohongCopynumFile(Xdir_1M + f, Xiaohong_segments)
    return Xiaohong_segments

def writeBalanceFiles(isegs, patient, all_samples):
    for sample in all_samples:
        outfile = open(outdir + sample + "_balanced_calls.txt", "w")
        outfile.write("Patient\tSample\tChr\tStart\tEnd\tCall\n")
        for chr in isegs:
            for iseg in isegs[chr]:
                outfile.write(patient)
                outfile.write("\t" + sample.split("_")[1])
                outfile.write("\t" + chr)
                outfile.write("\t" + str(iseg[0]))
                outfile.write("\t" + str(iseg[1]))
                matches = iseg[2]
                nonmatches = iseg[3]
                sample_in_matches = False
                sample_in_nonmatches = False
                for match in matches:
                    if sample in match:
                        sample_in_matches = True
                        break
                for nonmatch in nonmatches:
                    if sample in nonmatch:
                        sample_in_nonmatches = True
                        break
                if sample_in_matches:
                    outfile.write("\tUnbalanced")
                elif not sample_in_nonmatches:
                    outfile.write("\tShort")
                elif len(matches)>0:
                    outfile.write("\tBalanced")
                else:
                    outfile.write("\tUnknown")
                outfile.write("\n")

def writeSummary(isegs, patient, all_samples, all_analyses):
    if not onlyonepatient and isfile(outdir + patient + "_analysis_overview.txt"):
        print("Skipping patient", patient, ": analysis already exists.")
        return
    overview = {}
    overview_bases = {}
    calls = ["TP", "FP", "UP", "TN", "FN", "UN", "NC", "S"]
    for sample in all_samples:
        overview[sample] = {}
        overview_bases[sample] = {}
        for analysis in all_analyses:
            overview[sample][analysis] = {}
            overview_bases[sample][analysis] = {}
            for call in calls:
                overview[sample][analysis][call] = 0
                overview_bases[sample][analysis][call] = 0
    summary_out = open(outdir + patient + "_full_analysis.txt", "w")
    summary_out.write("Patient")
    summary_out.write("\tSample")
    summary_out.write("\tchr")
    summary_out.write("\tstart")
    summary_out.write("\tend")
    for analysis in all_analyses:
        summary_out.write("\t" + analysis)
    summary_out.write("\n")
    for sample in all_samples:
        for chr in isegs:
            for iseg in isegs[chr]:
                summary_out.write(patient)
                summary_out.write("\t" + sample.split("_")[1])
                summary_out.write("\t" + chr)
                summary_out.write("\t" + str(iseg[0]))
                summary_out.write("\t" + str(iseg[1]))
                for analysis in all_analyses:
                    summary_out.write("\t" + iseg[4][sample][analysis])
                    overview[sample][analysis][iseg[4][sample][analysis]] += 1
                    overview_bases[sample][analysis][iseg[4][sample][analysis]] += iseg[1] - iseg[0]
                summary_out.write("\n")
    summary_out.close()

    overview_out = open(outdir + patient + "_analysis_overview.txt", "w")
    overview_out.write("Patient")
    overview_out.write("\tSample")
    overview_out.write("\tAnalysis")
    for call in calls:
        overview_out.write("\t" + call + "_num")
    for call in calls:
        overview_out.write("\t" + call + "_len")
    overview_out.write("\n")
    for sample in all_samples:
        for analysis in all_analyses:
            overview_out.write(patient)
            overview_out.write("\t" + sample.split("_")[1])
            overview_out.write("\t" + analysis)
            for call in calls:
                overview_out.write("\t" + str(overview[sample][analysis][call]))
            for call in calls:
                overview_out.write("\t" + str(overview_bases[sample][analysis][call]))
            overview_out.write("\n")
    overview_out.close()

def getCanonicalAscatCallsFor(patient):
    if patient=="163":
        return [("24735", "900", "tetraploid"), ("23740", "900", "tetraploid"), ("23743", "800", "tetraploid"), ("23749", "200", "tetraploid")]
    if patient=="184":
        return [("24325", "2500", "diploid"), ("24325", "2500", "tetraploid"), ("24328", "2500", "diploid"), ("24328", "2500", "tetraploid"), ("24331", "2500", "diploid"), ("24331", "2500", "tetraploid"), ("24334", "2500", "diploid"), ("24334", "2500", "tetraploid"), ("24340", "2500", "diploid"), ("24340", "2500", "tetraploid"), ]
    if patient=="396":
        return [("24445", "100", "diploid"), ("24445", "100", "tetraploid"), ("24448", "100", "diploid"), ("24448", "100", "tetraploid"), ("24448", "1200", "diploid"), ("24448", "1200", "tetraploid"), ("24448", "300", "diploid"), ("24448", "300", "tetraploid"), ("24457", "1000", "diploid"), ("24457", "1000", "tetraploid"), ("24461", "100", "tetraploid")]
    if patient=="1047":
        return[("24271", "2500", "diploid"), ("24274", "2500", "diploid"), ("24277", "2500", "diploid"), ("24280", "2500", "diploid"), ("24271", "2500", "tetraploid"), ("24274", "2500", "tetraploid"), ("24277", "2500", "tetraploid"), ("24280", "2500", "tetraploid"), ]
    print("Unknown patient")
    return []

def readAscatSegmentationFor(patient, canon):
    (sample, gamma, ploidy) = canon
    ascfile = open(gamma_outdir + "/pASCAT_input_g" + gamma + "/" + ploidy + "/" + patient + "_fcn_ascat_segments.txt", "r")
    segments = {}
    for line in ascfile:
        if "Chr" in line:
            continue
        (__, id, chr, start, end, __, __, __, nA, nB) = line.split()
        if sample not in id:
            continue
        if chr not in segments:
            segments[chr] = []
        nA = int(nA)
        nB = int(nB)
        if nA > nB:
            temp = nA
            nA = nB
            nB = temp
        call = "wt"
        if nA == nB:
            if nA==0:
                call = "Double_d"
            elif nA==1:
                call = "wt"
            else:
                call = "Balanced_gain"
        else:
            if nA==0 and nB == 1:
                call = "Loss"
            elif nA==0:
                call = "CNLOH"
            else:
                call = "Gain"
        segments[chr].append((int(start), int(end), call))
    return segments

def getCallFor(start, end, seg):
    (tstart, tend, call) = seg
    if start >= tstart and end <= tend:
        return call
    if start < tstart and end > tstart:
        if end - tstart > tstart - start:
            return "Mostly_" + call
    if start < tend and end > tend:
        if end - tend < tend - start:
            return "Mostly_" + call
    return ""

def writeComparison(Xsegs, Asegs, patient, sample, gamma, ploidy, isegs):
    compareout = open(outdir + "compare_xiaohong_to_ascat_" + patient + "_" + sample + "_g" + gamma + "_" + ploidy + ".txt", "w")
    compareout.write("Patient")
    compareout.write("\tSample")
    compareout.write("\tChr")
    compareout.write("\tStart")
    compareout.write("\tEnd")
    compareout.write("\tXCall")
    compareout.write("\tACall")
    compareout.write("\tLength")
    compareout.write("\tSame?")
    compareout.write("\n")
    for chr in isegs:
        sumdiff = 0
        for (start, end) in isegs[chr]:
            xcall = ""
            acall = ""
            if chr in Xsegs[patient][sample]:
                for xseg in Xsegs[patient][sample][chr]:
                    xcall = getCallFor(start, end, xseg)
                    if xcall != "":
                        break
            for aseg in Asegs[chr]:
                acall = getCallFor(start, end, aseg)
                if acall != "":
                    break
            if xcall == "":
                xcall = "wt?"
            
            compareout.write(patient)
            compareout.write("\t" + sample)
            compareout.write("\t" + chr)
            compareout.write("\t" + str(start))
            compareout.write("\t" + str(end))
            compareout.write("\t" + xcall)
            compareout.write("\t" + acall)
            compareout.write("\t" + str(end-start))
            compareout.write("\t" + str(acall in xcall))
            compareout.write("\n")
            if acall not in xcall:
                sumdiff += end-start
        compareout.write(patient)
        compareout.write("\t" + sample)
        compareout.write("\t" + chr)
        compareout.write("\t" + "--")
        compareout.write("\t" + "--")
        compareout.write("\t" + "--")
        compareout.write("\t" + "--")
        compareout.write("\t" + str(sumdiff))
        compareout.write("\n")
        
            

Xiaohong_segments = readAllXiaohongSegmentation()

for patient in ["163", "184", "396", "1047"]:
    isegs = getIsegsFromCopynumberFileFor(patient)
    canonical_ascat = getCanonicalAscatCallsFor(patient)
    for canon in canonical_ascat:
        Ascat_segments = readAscatSegmentationFor(patient, canon)
        writeComparison(Xiaohong_segments, Ascat_segments, patient, canon[0], canon[1], canon[2], isegs)
        
