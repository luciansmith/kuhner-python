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

import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages/pymix")

import numpy
#import operator
import mixture

import lucianSNPLibrary as lsl

onlysomepatients = True
somepatients = ["521"]

VAFdir = "VAFclusters/"
#outdir = "VAFclusters_histograms/"
#outdir_low = "VAFclusters_histograms_low_oneplus_521/"

#if not path.isdir(outdir):
#    mkdir(outdir)
#if not path.isdir(outdir_low):
#    mkdir(outdir_low)

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
                #print("Switched directions: going down at", key, distance, localmax, val)
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
                #print("Switched directions: going up at", key, distance, localmax, val)
                direction="up"
                localmax = val
                maxkey = key
    return ret

(patientSampleMap, samplePatientMap) = lsl.getPatientSampleMap(dipvtet_file="calling_evidence_odds.tsv")
deletions, CNVs = lsl.loadDeletionsAndCNVs(samplePatientMap)

summary = open("summary_smoothed_and_fit.tsv", "w")
summary.write("Patient")
summary.write("\tSample")
summary.write("\tnPoints")
summary.write("\tCall")
summary.write("\tGroup")
#summary.write("\tMean x2")
#summary.write("\tStdev x2")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\n")

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
            if len(data[label][call]) < 100:
                #Skip groups with fewer than 100 VAFs.
                continue
            filename = patient + "_" + sample + "_" + makeFilename(label) + "_" + str(call[0]) + "_" + str(call[1]) + "_hist.png"
            hist = lsl.createPrintAndSaveHistogram(data[label][call], filename, 0.001, xdata="VAF", savefig=False, show=False)
            mean = numpy.mean(data[label][call])
            stdev = numpy.std(data[label][call])
            histmaxes = getHistMaxes(hist)
            print(patient, sample, label, call)
            ###THIS IS WHERE YOU FIND THE HISTOGRAM PEAKS###
            ##Data:  data[label][call]
            ##Peaks:  histmaxes
            ##Peak heights:  hist[histmaxes[n]]
            ##Stdev:  stdev

            emdata = mixture.DataSet()           
            emdata.fromList(data[label][call])
            numpeaks = len(histmaxes)
            gaussian_objects = []
            weights = []
            for i in xrange(numpeaks):
              n = mixture.NormalDistribution(histmaxes[i],stdev)
              gaussian_objects.append(n)
              weights.append(hist[histmaxes[i]])
            totweight = float(sum(weights))
            weights = [x/totweight for x in weights]
            mymix = mixture.MixtureModel(numpeaks,weights,gaussian_objects)
            # print "Before",mymix
            mymix.EM(emdata,40,0.1)
            # print "After",mymix
            print("Number of peaks=",mymix.G)
            for i in range(mymix.G):
              print(mymix.pi[i], mymix.components[i])

            summary.write(patient)
            summary.write("\t" + sample)
            summary.write("\t" + str(len(data[label][call])))
            summary.write("\t" + str(call))
            summary.write("\t" + label)
#            summary.write("\t" + str(mean*2))
#            summary.write("\t" + str(stdev*2))
            for i in range(mymix.G):
                histmax = histmaxes[i]
                summary.write("\t" + str(histmax*2))
#                whichsum.write("\t" + str(hist[histmax]))
                summary.write("\t" + str(2*mymix.components[i].distList[0].mu))
                summary.write("\t" + str(mymix.pi[i]*len(data[label][call])))
            summary.write("\n")

summary.close()
