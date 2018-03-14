#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:43:10 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

import lucianSNPLibrary as lsl

CNdir = "nonintegerCNs/"
balanced_dir = "balanced_calls/"

onlysomepatients = False
somepatients = ["997"]

firstpatients = ["17", "42", "55", "59", "74", "43", "184", "163", "396", "1047"]


def readBalancedCalls(patient, sample):
    balfile = open(balanced_dir + patient + "_" + sample + "_balanced_calls.tsv", "r")
    balcalls = {}
    for line in balfile:
        if "Chr" in line:
            continue
        (fpatient, fsample, chr, start, end, call) = line.rstrip().split()
        if chr not in balcalls:
            balcalls[chr] = {}
        balcalls[chr][(int(start), int(end))] = call
    return balcalls


def readAmbiguousCallsFromASCAT(f):
    cnfile = open(CNdir + f, "r")
    unbal_calls = {}
    for line in cnfile:
        if "patient" in line:
            continue
        (patient, sample, chr, start, end, rawA, rawB, intA, intB) = line.split()
        if intA == "NA" or intB == "NA":
            continue
        rawA = float(rawA)
        rawB = float(rawB)
        if abs(rawA-rawB) < 0.000001 and intA != intB:
            if chr not in unbal_calls:
                unbal_calls[chr] = []
            unbal_calls[chr].append((int(start), int(end), rawA))
    return unbal_calls

def compareAndReport(unbal_calls, bal_calls, patient, sample, ploidy, allbal, allunbal):
    for chr in unbal_calls:
        for ub_seg in unbal_calls[chr]:
            for b_seg in bal_calls[chr]:
                if bal_calls[chr][b_seg] == "Unknown":
                    continue
                (bstart, bend) = b_seg
                (ubstart, ubend, rawA) = ub_seg
                if bstart >= ubstart and bend <= ubend:
                    if bal_calls[chr][b_seg] == "Balanced":
                        allbal.append(rawA)
#                        print("Balanced:", chr, ub_seg)
                    else:
                        allunbal.append(rawA)
#                        print("Unbalanced:", chr, ub_seg)

allbal = []
allunbal = []
files = []
for (__, __, f) in walk(CNdir):
    files += f
for f in files:
    if "nonint" not in f:
        continue
    (patient, sample, ploidy) = f.split("_")[0:3]
    if onlysomepatients and patient not in somepatients:
        continue
    unbal_calls = readAmbiguousCallsFromASCAT(f)
    bal_calls = readBalancedCalls(patient, sample)
    print ("Comparing", patient, sample, ploidy)
    compareAndReport(unbal_calls, bal_calls, patient, sample, ploidy, allbal, allunbal)

print("All balanced values:")
lsl.createPrintAndSaveHistogram(allbal, "", 0.001)
print("All unbalanced values:")
lsl.createPrintAndSaveHistogram(allunbal, "", 0.001)
