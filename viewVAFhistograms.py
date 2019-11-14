#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 12:45:13 2019

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

import numpy
import operator

import lucianSNPLibrary as lsl

onlysomepatients = True
somepatients = ["521"]

VAFdir = "VAFclusters_kanika/"
outdir = "VAFclusters_histograms_kanika/"
outdir_low = "VAFclusters_histograms_low_kanika/"

if not path.isdir(outdir):
    mkdir(outdir)
if not path.isdir(outdir_low):
    mkdir(outdir_low)

VAFfiles = []
for __, _, files in walk(VAFdir):
    VAFfiles += files

def sortLabels(labels):
    newlist = []
    sublist = []
    for label in labels:
        if len(label.split(", '"))==1:
            sublist.append(label)
    sublist.sort()
    newlist.extend(sublist)
    for n in range(9,1,-1):
#    for n in range(2,9):
        sublist = []
        for label in labels:
            if len(label.split(", '"))==n:
                sublist.append(label)
        sublist.sort()
        newlist.extend(sublist)
    return newlist

def getCNVCall(patient, sample, chrom, pos, CNVs):
    if patient not in CNVs:
        assert(False)
        return (-1, -1)
    if sample not in CNVs[patient]:
        assert(False)
        return (-1, -1)
    if chrom not in CNVs[patient][sample]:
        print("No chromosome", str(chrom), "found.")
        assert(False)
        return (-1, -1)
    for (start, end, call) in CNVs[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return call
    return (-1, -1)

def makeFilename(label):
    elements = eval(label)
    ret = ""
    for element in elements:
        ret += element + "-"
    return ret[:-1]

def getHistMaxes(hist):
    ret = []
    keylist = list(hist.keys())
    keylist.sort()
    localmax = 0
    maxkey = keylist[0]
    localmin = 0
    direction = "up"
    distance = 0
    for key in keylist:
        val = hist[key]
        if direction=="up":
            if val > localmax:
                localmax = val
                maxkey = key
                distance = 0
            else:
                distance += 1
            if distance >= 40 and localmax-val > 0.25:
                print("Switched directions: going down at", key, distance, localmax, val)
                direction = "down"
                ret.append(maxkey)
                localmin = val
        elif direction=="down":
            if val < localmin:
                localmin = val
                distance = 0
            else:
                distance += 1
            if distance >= 40 and val - localmin > 0.25:
                print("Switched directions: going up at", key, distance, localmax, val)
                direction="up"
                localmax = val
                maxkey = key
    return ret

(patientSampleMap, samplePatientMap) = lsl.getPatientSampleMap(dipvtet_file="calling_evidence_odds.tsv")
deletions, CNVs = lsl.loadDeletionsAndCNVs(samplePatientMap)


summary = open(outdir + "summary.tsv", "w")
summary.write("Patient")
summary.write("\tSample")
summary.write("\tnPoints")
summary.write("\tCall")
summary.write("\tGroup")
summary.write("\tMean x2")
summary.write("\tStdev x2")
summary.write("\tHistMax x2")
summary.write("\tHistMax height")
summary.write("\tHistMax x2")
summary.write("\tHistMax height")
summary.write("\tHistMax x2")
summary.write("\tHistMax height")
summary.write("\tHistMax x2")
summary.write("\tHistMax height")
summary.write("\tHistMax x2")
summary.write("\tHistMax height")
summary.write("\tHistMax x2")
summary.write("\tHistMax height")
summary.write("\n")

summary_low = open(outdir_low + "summary.tsv", "w")
summary_low.write("Patient")
summary_low.write("\tSample")
summary_low.write("\tnPoints")
summary_low.write("\tCall")
summary_low.write("\tGroup")
summary_low.write("\tMean x2")
summary_low.write("\tStdev x2")
summary_low.write("\tHistMax x2")
summary_low.write("\tHistMax height")
summary_low.write("\tHistMax x2")
summary_low.write("\tHistMax height")
summary_low.write("\tHistMax x2")
summary_low.write("\tHistMax height")
summary_low.write("\tHistMax x2")
summary_low.write("\tHistMax height")
summary_low.write("\tHistMax x2")
summary_low.write("\tHistMax height")
summary_low.write("\tHistMax x2")
summary_low.write("\tHistMax height")
summary_low.write("\n")

for file in VAFfiles:
    if "_VAFs" not in file:
        continue
    (patient, sample) = file.split("_")[0:2]
    if onlysomepatients and patient not in somepatients:
        continue
    data = {}
    for line in open(VAFdir + file, "r"):
        lvec = line.rstrip().split("\t")
        if "Patient" in line:
            labels = lvec[5:]
            for group in labels:
                data[group] = {}
            continue
        if lvec[2]=="23" or lvec[2]=="24":
            continue
        call = getCNVCall(patient, sample, lvec[2], int(lvec[3]), CNVs)
        if call==(-1, -1):
            continue
        for n in range(5, len(lvec)):
            group = labels[n-5]
            if lvec[n] != "":
                if call not in data[group]:
                    data[group][call] = []
                data[group][call].append(float(lvec[n]))
    for label in data:
        for call in data[label]:
            if len(data[label][call]) == 0:
                continue
            direc = outdir
            whichsum = summary
            if len(data[label][call]) < 100:
                direc = outdir_low
                whichsum = summary_low
            filename = direc + patient + "_" + sample + "_" + makeFilename(label) + "_" + str(call[0]) + "_" + str(call[1]) + "_hist.png"
            hist = lsl.createPrintAndSaveHistogram(data[label][call], filename, 0.001, xdata="VAF", savefig=True, show=False)
            mean = numpy.mean(data[label][call])
            stdev = numpy.std(data[label][call])
            histmaxes = getHistMaxes(hist)
            whichsum.write(patient)
            whichsum.write("\t" + sample)
            whichsum.write("\t" + str(len(data[label][call])))
            whichsum.write("\t" + str(call))
            whichsum.write("\t" + label)
            whichsum.write("\t" + str(mean*2))
            whichsum.write("\t" + str(stdev*2))
            for histmax in histmaxes:
                whichsum.write("\t" + str(histmax*2))
                whichsum.write("\t" + str(hist[histmax]))
            whichsum.write("\n")
            print(patient, sample, label, call, direc)
            
summary.close()
summary_low.close()

