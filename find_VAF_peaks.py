#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 13:09:24 2018

@author: lpsmith
"""


from __future__ import division
from os import walk
from os import path
from os import mkdir
from shutil import copytree

VAF_2v4dir = "VAF_2v4_hist/"

onlysomepatients = False
somepatients = ["997"]

VAFfile = open (VAF_2v4dir + "2v4_peaks.tsv", "w")
VAFfile.write("Patient")
VAFfile.write("\tSample")
VAFfile.write("\tVAF")
VAFfile.write("\tProbability")
VAFfile.write("\tPeak or Valley")
VAFfile.write("\tnVAFs")
VAFfile.write("\n")


infiles = []
for (_, _, f) in walk(VAF_2v4dir):
    infiles += f
for file in infiles:
    if "VAF_hist.tsv" not in file:
        continue
    (patient, sample) = file.split("_")[0:2]
    if onlysomepatients and patient not in somepatients:
        continue
    hist = []
    nVAFs = 0
    for line in open(VAF_2v4dir + file, "r"):
        if "VAF" in line:
            nVAFs = line.split()[3]
            continue
        (VAF, prob) = line.split()
        VAF = float(VAF)
        prob = float(prob)
        hist.append((VAF, prob))
    hist.sort()
    peakindex = 0
    peakprob = 0
    allpeaks = []
    slope = "pos"
    for index in range(0,len(hist)):
        (VAF, prob) = hist[index]
        if prob > peakprob:
            if slope=="neg":
                allpeaks.append((hist[peakindex][0], hist[peakindex][1], "Valley"))
            peakindex = index
            peakprob = prob
            slope = "pos"
        if peakindex + 20 < index:
            if slope=="pos":
                allpeaks.append((hist[peakindex][0], hist[peakindex][1], "Peak"))
            slope = "neg"
            peakindex += 1
            peakprob = hist[peakindex][1]
    for peak in allpeaks:
        VAFfile.write(patient)
        VAFfile.write("\t" + sample)
        VAFfile.write("\t" + str(peak[0]))
        VAFfile.write("\t" + str(peak[1]))
        VAFfile.write("\t" + str(peak[2]))
        VAFfile.write("\t" + nVAFs)
        VAFfile.write("\n")

VAFfile.close()

