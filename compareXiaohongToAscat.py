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
from os import mkdir
from shutil import copy2 as copy
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

processed_ascat_dir = "noninteger_processed_CNs/"
gamma_outdir = "gamma_test_output/"
balanced_dir = "balanced_calls/"
patient_dir = gamma_outdir + "pASCAT_input_g500/"
outdir = "Xiaohong_pASCAT_compare/"
jamboree_dir = "jamboree_files/"
gamma_list = ["100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500"]

use500 = True
onlysomepatients = False
#somepatients = ["55", "59", "74", "591", "595", "597"]
#somepatients = ["403", "512", "568", "852", "572"]
somepatients = ["772"]

if not path.isdir(outdir):
    mkdir(outdir)
if not path.isdir(jamboree_dir):
    mkdir(jamboree_dir)

def getIsegsFromCopynumberFileFor(patient):
    #This function reads the lowest-gamma-value copynumber file it can, reads it, and returns it as 'isegs'.
    isegs = {}
    #glist = ("100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000")
    glist = ("500",)
    if patient=="772":
        glist = ("550",)
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
#        if nbaf <10:
#            continue
        start = int(start)
        end = int(end)
        if chr not in isegs:
            isegs[chr] = []
        isegs[chr].append([start, end, {}])
    return isegs

def readBalancedUnbalancedAndStoreInIsegs(isegs, patient):
    files = []
    for (__, __, f) in walk(balanced_dir):
        files += f
    for f in files:
        if "balanced_calls" not in f:
            continue
        (ipatient, sample, __, __) = f.split("_")
        if patient != ipatient:
            continue
        balfile = open(balanced_dir + f, "r")
        print("Reading file", f)
        for line in balfile:
            if "Chr" in line:
                continue
            (fpatient, fsample, chr, start, end, call) = line.rstrip().split()
            if fpatient != patient:
                print("Wrong patient!", f, patient, fpatient)
                continue
            start = int(start)
            end = int(end)
            if chr in isegs:
                for iseg in isegs[chr]:
                    if iseg[0] == start and iseg[1] == end:
                        iseg[2][sample] = call
                        break


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

def addCentromereToXiaohong(xsegs):
    for patient in xsegs:
        for sample in xsegs[patient]:
            for chr in xsegs[patient][sample]:
                (start, end) = lsl.getChromosomeCentromere(int(chr))
                xsegs[patient][sample][chr].append([start, end, "Centromere"])

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
    addCentromereToXiaohong(Xiaohong_segments)
    return Xiaohong_segments

def readAscatSegmentationFor(patient, canon):
    (sample, gamma, ploidy) = canon
    ascfile = open(processed_ascat_dir + patient + "_" + sample + "_g" + gamma + "_" + ploidy + "_nonint_CNs.txt", "r")
    segments = {}
    for line in ascfile:
        if "atient" in line:
            continue
        (id1, id2, chr, start, end, __, __, nA, nB) = line.split()
        if sample not in id2:
            continue
        if chr not in segments:
            segments[chr] = []
        if nA=="NA" or nB=="NA":
            call = "Unknown"
        else:
            nA = int(nA)
            nB = int(nB)
            nsum = nA + nB
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
                if nA==0 or nB == 0:
                    if nsum == 1:
                        call = "Loss"
                    elif nsum == 2:
                        call = "CNLOH"
                    else:
                        call = "LOH_Gain"
                else:
                    call = "Gain"
        segments[chr].append((int(start), int(end), call))
    return segments

def getCallFor(start, end, seg):
    (tstart, tend, call) = seg
    overlap = max(0, (min(end, tend) - max(start, tstart) + 1)) / (end-start+1)
    return (call, overlap)

def addWtToCalls(calls):
    totcalls = 0
    for call in calls:
        totcalls += call[1]
    if totcalls < 0.98:
        calls.append(("wt?", 1-totcalls))

def addUnknownToCalls(calls):
    totcalls = 0
    for call in calls:
        totcalls += call[1]
    if totcalls < 0.98:
        #assert(False) #When only using one gamma, ASCAT calls should always overlap perfectly
        calls.append(("Unknown", 1-totcalls))


def getSummary(calls):
    if len(calls)==0:
        assert(False)
        #if wtNotCalled:
        #    return "wt?"
        return "??"
    totcalls = 0
    summary = {}
    for call in calls:
        totcalls += call[1]
        callname = call[0]
        if callname not in summary:
            summary[callname] = 0
        summary[callname] += call[1]
    if len(summary.keys()) == 1:
        return list(summary.keys())[0] #'keys' is an iterator: cast to a list and take the first (and only) entry
    ret = ""
    for sumtype in summary:
        ret += sumtype + "_" + str(round(summary[sumtype], 2)) + "_"
    return ret

def getMatches(bcall, calls):
    valid = 0
    contradicted = 0
    for call in calls:
        if bcall == "Balanced":
            if "wt" in call[0] or "Double_d" in call[0] or "Balanced_gain" in call[0] or "Centromere" in call[0]:
                valid += call[1]
            else:
                contradicted += call[1]
        elif bcall == "Unbalanced":
            if "Gain" in call[0] or "Loss" in call[0] or "LOH" in call[0] or "Centromere" in call[0]:
                valid += call[1]
            else:
                contradicted += call[1]
    if valid == 0 and contradicted == 0:
        return "Unknown"
    if valid > contradicted:
        return "Validated"
    return "Contradicted"


def getMatchTotal(xcalls, acalls):
    matchtotal = 0
    for acall in acalls:
        if acall[0] == "Unknown":
            matchtotal += acall[1]
            continue
        for xcall in xcalls:
            if acall[0] in xcall[0]:
                matchtotal += min(xcall[1], acall[1])
            elif acall[0] == "LOH_Gain" and ("Gain" in xcall[0] or "LOH" in xcall[0]):
                matchtotal += min(xcall[1], acall[1])
            elif xcall[0] == "Centromere":
                matchtotal += min(xcall[1], acall[1])
#    if len(acalls) > 1:
#        print(acalls, xcalls, matchtotal)
    if matchtotal > 0.95:
        return 1
    return matchtotal

def writeComparison(Xsegs, Asegs, patient, sample, gamma, ploidy, isegs):
    compareout = open(outdir + patient + "_" + sample + "_g" + gamma + "_" + ploidy + "_xiaohong_to_ascat_compare.tsv", "w")
    compareout.write("Patient")
    compareout.write("\tSample")
    compareout.write("\tChr")
    compareout.write("\tStart")
    compareout.write("\tEnd")
    compareout.write("\tXCall")
    compareout.write("\tACall")
    compareout.write("\tLength")
    compareout.write("\tLength mismatch")
    compareout.write("\tBalanced")
    compareout.write("\tXCall_validated")
    compareout.write("\tACall_validated")
    compareout.write("\n")
    for chr in isegs:
        sumdiff = 0
        for (start, end, balanced_calls) in isegs[chr]:
            xcalls = []
            acalls = []
            if sample in Xsegs[patient] and chr in Xsegs[patient][sample]:
                for xseg in Xsegs[patient][sample][chr]:
                    xcall = getCallFor(start, end, xseg)
                    if xcall[1] > 0:
                        xcalls.append(xcall)
            for aseg in Asegs[chr]:
                acall = getCallFor(start, end, aseg)
                if acall[1] > 0:
                    acalls.append(acall)
            addWtToCalls(xcalls)
            addUnknownToCalls(acalls)
            xcall = getSummary(xcalls)
            acall = getSummary(acalls)
            
            x_matches_balanced = getMatches(balanced_calls[sample], xcalls)
            a_matches_balanced = getMatches(balanced_calls[sample], acalls)
            
            length = end-start+1
            a_x_mismatch = int(length - (length*getMatchTotal(xcalls, acalls)))

            compareout.write(patient)
            compareout.write("\t" + sample)
            compareout.write("\t" + chr)
            compareout.write("\t" + str(start))
            compareout.write("\t" + str(end))
            compareout.write("\t" + xcall)
            compareout.write("\t" + acall)
            compareout.write("\t" + str(length))
            compareout.write("\t" + str(a_x_mismatch))
            compareout.write("\t" + balanced_calls[sample])
            compareout.write("\t" + x_matches_balanced)
            compareout.write("\t" + a_matches_balanced)
            compareout.write("\n")
            sumdiff += a_x_mismatch
        compareout.write(patient)
        compareout.write("\t" + sample)
        compareout.write("\t" + chr)
        compareout.write("\t" + "--")
        compareout.write("\t" + "--")
        compareout.write("\t" + "--")
        compareout.write("\t" + "--")
        compareout.write("\t" + str(sumdiff))
        compareout.write("\n")

def copyDataForJamboree(patient, data, gamma, ploidy):
    pseg_dir = "pASCAT_segmentation/"
    if not path.isdir(jamboree_dir + pseg_dir):
        mkdir(jamboree_dir + pseg_dir)
    if not path.isdir(jamboree_dir + balanced_dir):
        mkdir(jamboree_dir + balanced_dir)
    xcompare_dir = "xiaohong_compare/"
    if not path.isdir(jamboree_dir + xcompare_dir):
        mkdir(jamboree_dir + xcompare_dir)
    infile = open(gamma_outdir + "pASCAT_input_g" + gamma + "/" + ploidy  + "/" + patient + "_fcn_ascat_segments.txt", "r")
    outfile = open(jamboree_dir + pseg_dir + patient + "_" + sample + "_g" + gamma + "_" + ploidy + "_segments.tsv", "w")
    for line in infile:
        if "SampleID" in line or patient + "_" + sample in line:
            outfile.write(line)
    copy(outdir + patient + "_" + sample + "_g" + gamma + "_" + ploidy + "_xiaohong_to_ascat_compare.tsv", jamboree_dir + xcompare_dir)
    copy(balanced_dir + patient + "_" + sample + "_balanced_calls.tsv", jamboree_dir + balanced_dir)
    

Xiaohong_segments = readAllXiaohongSegmentation()

files = []
for (__, __, f) in walk(patient_dir):
    files += f
for f in files:
    if "copynumber_segments" not in f:
        continue
    patient = f.split("_")[0]
    if onlysomepatients and patient not in somepatients:
        continue
#for patient in ["43"]:
    isegs = getIsegsFromCopynumberFileFor(patient)
    readBalancedUnbalancedAndStoreInIsegs(isegs, patient)
    if use500:
        canonical_ascats = lsl.getSingleGammaCallsFor(patient, "500")
    else:
        canonical_ascats = lsl.getCanonicalAscatCallsFor(patient)
    for canon in canonical_ascats:
        sample = canon[0]
        gamma = canon[1]
        ploidy = canon[2]
        if gamma == "None":
            continue
        Ascat_segments = readAscatSegmentationFor(patient, canon)
        writeComparison(Xiaohong_segments, Ascat_segments, patient, sample, gamma, ploidy, isegs)
        if "N" not in sample and int(sample) < 23341 and sample != "19578":
            continue
        copyDataForJamboree(patient, sample, gamma, ploidy)

