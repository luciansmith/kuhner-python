# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:24:27 2018

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

onlysomepatients = False
somepatients = ["74"]


xiaocompare_dir = "Xiaohong_pASCAT_compare/"

sameness = ["same", "different"]
categories = ["Contradicted", "Unknown", "Validated"]

summaries = {}

def addtoSummaries(summaries, patient, sample, gamma, ploidy, xcall, acall, length, issame, xvalid, avalid):
    length = int(length)
    if patient not in summaries:
        summaries[patient] = {}
    if sample not in summaries[patient]:
        summaries[patient][sample] = {}
        summaries[patient][sample]["gploidy"] = (gamma, ploidy)
        for same in sameness:
            summaries[patient][sample][same] = {}
            if same=="same":
                for category in categories:
                    summaries[patient][sample][same][category] = 0
            else:
                summaries[patient][sample][same]["Unknown"] = 0
                for xv in ["Contradicted", "Validated"]:
                    summaries[patient][sample][same][xv] = {}
                    for av in ["Contradicted", "Validated"]:
                        summaries[patient][sample][same][xv][av] = 0
    if issame == "True":
        same = "same"
        summaries[patient][sample][same][xvalid] += length
    elif issame == "False":
        same = "different"
        if xvalid=="Unknown":
            summaries[patient][sample][same][xvalid] += length
        else:
            summaries[patient][sample][same][xvalid][avalid] += length
    else:
        print("Unknown 'is same':", issame)

def writeHeader(sumout):
    sumout.write("Patient")
    sumout.write("\tSample")
    sumout.write("\tGamma")
    sumout.write("\tPloidy")
    sumout.write("\tSame, validated")
    sumout.write("\tSame, contradicted")
    sumout.write("\tSame, unknown")
    sumout.write("\tSame total")
    sumout.write("\tDifferent, both valid")
    sumout.write("\tDifferent, both contradicted")
    sumout.write("\tDifferent, unknown")
    sumout.write("\tDifferent, X valid; A contradicted")
    sumout.write("\tDifferent, A valid; X contradicted")
    sumout.write("\tDifferent total")
    sumout.write("\n")

def writeLine(sumout, psummaries, patient, sample):
    for sample in psummaries:
        sumout.write(patient)
        sumout.write("\t" + sample)
        sumout.write("\t" + str(psummaries[sample]["gploidy"][0]))
        sumout.write("\t" + str(psummaries[sample]["gploidy"][1]))
        sumout.write("\t" + str(psummaries[sample]["same"]["Validated"]))
        sumout.write("\t" + str(psummaries[sample]["same"]["Contradicted"]))
        sumout.write("\t" + str(psummaries[sample]["same"]["Unknown"]))
        total = 0
        for category in psummaries[sample]["same"]:
            total += psummaries[sample]["same"][category]
        sumout.write("\t" + str(total))
        sumout.write("\t" + str(psummaries[sample]["different"]["Validated"]["Validated"]))
        sumout.write("\t" + str(psummaries[sample]["different"]["Contradicted"]["Contradicted"]))
        sumout.write("\t" + str(psummaries[sample]["different"]["Unknown"]))
        sumout.write("\t" + str(psummaries[sample]["different"]["Validated"]["Contradicted"]))
        sumout.write("\t" + str(psummaries[sample]["different"]["Contradicted"]["Validated"]))
        total = 0
        for category in psummaries[sample]["different"]:
            if category == "Unknown":
                total += psummaries[sample]["different"][category]
            else:
                for cat2 in psummaries[sample]["different"][category]:
                    total += psummaries[sample]["different"][category][cat2]
        sumout.write("\t" + str(total))
        sumout.write("\n")


def saveSummaries(summaries):
    sumout = open(xiaocompare_dir + "xiaocompare_summary.tsv", "w")
    jsumout = open(xiaocompare_dir + "xiaocompare_jonly_summary.tsv", "w")
    writeHeader(sumout)
    writeHeader(jsumout)
    for patient in summaries:
        for sample in summaries[patient]:
            writeLine(sumout, summaries[patient], patient, sample)
            if int(sample) >= 23341:
                writeLine(jsumout, summaries[patient], patient, sample)
    sumout.close()
    jsumout.close()
        
        

files = []
for (__, __, f) in walk(xiaocompare_dir):
    files += f
for f in files:
    if "xiaohong_to_ascat_compare" not in f:
        continue
    (patient, sample, gamma, ploidy) = f.split("_")[0:4]
    if onlysomepatients and patient not in somepatients:
        continue
    compare = open(xiaocompare_dir + f, "r")
    for line in compare:
        if "Patient" in line:
            continue
        lvec = line.split()
        if len(lvec) == 8:
            continue #Recap-chromsome line
        (patient, sample, chr, start, end, xcall, acall, length, same, balanced, xvalid, avalid) = lvec
        addtoSummaries(summaries, patient, sample, gamma, ploidy, xcall, acall, length, same, xvalid, avalid)

saveSummaries(summaries)
    
