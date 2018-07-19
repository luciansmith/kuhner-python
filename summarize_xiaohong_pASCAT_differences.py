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
ascat_dir = "gamma_test_output/"
best_dir = "best_analyses/"

sameness = ["same", "different"]
categories = ["Contradicted", "Unknown", "Validated"]

summaries = {}

def addtoSummaries(summaries, patient, sample, gamma, ploidy, xcall, acall, length, difflen, xvalid, avalid):
    length = int(length)
    if patient not in summaries:
        summaries[patient] = {}
    samGamPloidy = (sample, gamma, ploidy)
    if samGamPloidy not in summaries[patient]:
        summaries[patient][samGamPloidy] = {}
        summaries[patient][samGamPloidy]["gploidy"] = (gamma, ploidy)
        for same in sameness:
            summaries[patient][samGamPloidy][same] = {}
            if same=="same":
                for category in categories:
                    summaries[patient][samGamPloidy][same][category] = 0
            else:
                summaries[patient][samGamPloidy][same]["Unknown"] = 0
                for xv in ["Contradicted", "Validated"]:
                    summaries[patient][samGamPloidy][same][xv] = {}
                    for av in ["Contradicted", "Validated"]:
                        summaries[patient][samGamPloidy][same][xv][av] = 0
    difflen = int(difflen)
    same = "same"
    summaries[patient][samGamPloidy][same][xvalid] += length - difflen
    same = "different"
    if xvalid=="Unknown":
        summaries[patient][samGamPloidy][same][xvalid] += difflen
    else:
        summaries[patient][samGamPloidy][same][xvalid][avalid] += difflen

def writeHeader(sumout):
    sumout.write("Patient")
    sumout.write("\tSample")
    sumout.write("\tGamma")
    sumout.write("\tPloidy")
    sumout.write("\tPloidyVal")
    sumout.write("\tPurity")
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
    sumout.write("\tX accuracy")
    sumout.write("\tA accuracy")
    sumout.write("\n")

def writeLine(sumout, psummaries, patient, samGamPloidy, ploidyval, purityval, accuracies):
    sumout.write(patient)
    sumout.write("\t" + samGamPloidy[0])
    sumout.write("\t" + samGamPloidy[1])
    sumout.write("\t" + samGamPloidy[2])
    sumout.write("\t" + ploidyval)
    sumout.write("\t" + purityval)
    sumout.write("\t" + str(psummaries[samGamPloidy]["same"]["Validated"]/1000000))
    sumout.write("\t" + str(psummaries[samGamPloidy]["same"]["Contradicted"]/1000000))
    sumout.write("\t" + str(psummaries[samGamPloidy]["same"]["Unknown"]/1000000))
    total = 0
    for category in psummaries[samGamPloidy]["same"]:
        total += psummaries[samGamPloidy]["same"][category]
    sumout.write("\t" + str(total/1000000))
    sumout.write("\t" + str(psummaries[samGamPloidy]["different"]["Validated"]["Validated"]/1000000))
    sumout.write("\t" + str(psummaries[samGamPloidy]["different"]["Contradicted"]["Contradicted"]/1000000))
    sumout.write("\t" + str(psummaries[samGamPloidy]["different"]["Unknown"]/1000000))
    sumout.write("\t" + str(psummaries[samGamPloidy]["different"]["Validated"]["Contradicted"]/1000000))
    sumout.write("\t" + str(psummaries[samGamPloidy]["different"]["Contradicted"]["Validated"]/1000000))
    total = 0
    for category in psummaries[samGamPloidy]["different"]:
        if category == "Unknown":
            total += psummaries[samGamPloidy]["different"][category]
        else:
            for cat2 in psummaries[samGamPloidy]["different"][category]:
                total += psummaries[samGamPloidy]["different"][category][cat2]
    sumout.write("\t" + str(total/1000000))
    sumout.write("\t" + accuracies[0])
    sumout.write("\t" + accuracies[1])
    sumout.write("\n")

def getPurityAndPloidyVal(patient, samGamPloidy):
    purityval = "??"
    ploidyval = "??"
    (sample, gamma, ploidy) = samGamPloidy
    purityfile = open(ascat_dir + "pASCAT_input_" + gamma + "/" + ploidy + "/" + patient + "_fcn_ascat_cont.txt", "r")
    for line in purityfile:
        if patient + "_" + sample in line:
            purityval = line.split()[1]
    purityfile.close()
    ploidyfile = open(ascat_dir + "pASCAT_input_" + gamma + "/" + ploidy + "/" + patient + "_fcn_ascat_ploidy.txt", "r")
    for line in ploidyfile:
        if patient + "_" + sample in line:
            ploidyval = line.split()[1]
    ploidyfile.close()
    return (purityval, ploidyval)

def getAccuracies(patient, samGamPloidy):
    xaccuracy = "??"
    aaccuracy = "??"
    (sample, gamma, ploidy) = samGamPloidy
    bestfile = open(best_dir + patient + "_best.tsv", "r")
    for line in bestfile:
        (bpatient, bsample, bconstraint, bcompare, baccuracy, bgamma) = line.split()[0:6]
        if bsample != sample:
            continue
        if bcompare != "by_length":
            continue
        if bconstraint == "Xiaohong":
            xaccuracy = baccuracy
        elif bconstraint == ploidy:
            aaccuracy = baccuracy
    bestfile.close()
    return (xaccuracy, aaccuracy)

def saveSummaries(summaries):
    sumout = open(xiaocompare_dir + "xiaocompare_summary.tsv", "w")
    jsumout = open(xiaocompare_dir + "xiaocompare_jonly_summary.tsv", "w")
    writeHeader(sumout)
    writeHeader(jsumout)
    for patient in summaries:
        for samGamPloidy in summaries[patient]:
            (purityval, ploidyval) = getPurityAndPloidyVal(patient, samGamPloidy)
            accuracies = getAccuracies(patient, samGamPloidy)
            writeLine(sumout, summaries[patient], patient, samGamPloidy, ploidyval, purityval, accuracies)
            if "N" in samGamPloidy[0] or int(samGamPloidy[0]) >= 23341 or samGamPloidy[0] == 19578:
                writeLine(jsumout, summaries[patient], patient, samGamPloidy, ploidyval, purityval, accuracies)
    sumout.close()
    jsumout.close()
    copy(xiaocompare_dir + "xiaocompare_jonly_summary.tsv", "jamboree_files/xiaohong_compare/")
        
        

files = []
for (__, __, f) in walk(xiaocompare_dir):
    files += f
for f in files:
    if "xiaohong_to_ascat_compare" not in f:
        continue
    (patient, sample, gamma, ploidy) = f.split("_")[0:4]
    print("Summarizing", patient, sample, gamma, ploidy)
    if onlysomepatients and patient not in somepatients:
        continue
    compare = open(xiaocompare_dir + f, "r")
    for line in compare:
        if "Patient" in line:
            continue
        lvec = line.split()
        if len(lvec) == 8:
            continue #Recap-chromsome line
        (patient, sample, chr, start, end, xcall, acall, length, difflen, balanced, xvalid, avalid) = lvec
        addtoSummaries(summaries, patient, sample, gamma, ploidy, xcall, acall, length, difflen, xvalid, avalid)

saveSummaries(summaries)
    
