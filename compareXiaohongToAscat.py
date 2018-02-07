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

gamma_outdir = "gamma_test_output/"
outdir = "gamma_test_output/analysis_compare/"
gamma_list = ["0", "50", "100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000"]

onlyonepatient = False
onepatient = "163"

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
    overlap = max(0, (min(end, tend) - max(start, tstart))) / (end-start)
    return (call, overlap)

def getSummary(calls, wtNotCalled):
    if len(calls)==1:
        return calls[0][0]
    if len(calls)==0:
        if wtNotCalled:
            return "wt?"
        return "??"
    totcalls = 0
    summary = {}
    for call in calls:
        totcalls += call[1]
        callname = call[0]
        if callname not in summary:
            summary[callname] = 0
        summary[callname] += call[1]
    if totcalls < 0.8 and wtNotCalled:
        summary["wt"] = 1-totcalls
    if len(summary.keys()) == 1:
        return list(summary.keys())[0]
    ret = ""
    for sumtype in summary:
        ret += sumtype + "_" + str(round(summary[sumtype], 2)) + "_"
    return ret

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
            xcalls = []
            acalls = []
            if chr in Xsegs[patient][sample]:
                for xseg in Xsegs[patient][sample][chr]:
                    xcall = getCallFor(start, end, xseg)
                    if xcall[1] > 0:
                        xcalls.append(xcall)
            for aseg in Asegs[chr]:
                acall = getCallFor(start, end, aseg)
                if acall[1] > 0:
                    acalls.append(acall)
            xcall = getSummary(xcalls, True)
            acall = getSummary(acalls, False)
            
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
        