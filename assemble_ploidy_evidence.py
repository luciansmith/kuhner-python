#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 14:09:13 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os.path import isfile
from os import mkdir

#Process various input files to create an evidence file

initcall = "initial_calling_evidence.txt"
flow_summary = "flow_summary.txt"
call_summary = "Xiaohong_pASCAT_compare/xiaocompare_summary.tsv"
outfile = "calling_evidence.txt"

def writeNewLine(f, line, additions):
    f.write(line.rstrip())
    for addition in additions:
        f.write("\t" + str(addition))
    f.write("\n")


def readFlowSummary(flow):
    flowdata = {}
    for line in open(flow, "r"):
        if "Patient" in line:
            continue
        lvec = line.split()
        (progression, patient, ratio) = lvec[0:3]
        aneuploid_strs = lvec[3:len(lvec)]
        aneuploids = set()
        for aneuploid in aneuploid_strs:
            aneuploids.add(float(aneuploid))
        ratio = ratio.replace(':','::')
        flowdata[patient] = (ratio, progression, aneuploids)
    return flowdata

def readCallSummary(call):
    calldata = {}
    for line in open(call, "r"):
        if "Patient" in line:
            continue
        (patient, sample, gamma, ploidy, ploidyval, purval, __, __, __, __, __, __, __, __, __, __, Xacc, Aacc) = line.split()
        ploidyval = float(ploidyval)
        purval = float(purval)
        Aacc = float(Aacc)
        if patient not in calldata:
            calldata[patient] = {}
        if sample not in calldata[patient]:
            calldata[patient][sample] = {}
        calldata[patient][sample][ploidy] = (ploidyval, purval, Aacc)
    return calldata

def getDtMatches(patient, sample, flowdata, calldata):
    dmatch = "None"
    tmatch = "None"
    if patient in flowdata and patient in calldata:
        flowlist = flowdata[patient][2]
        if sample in calldata[patient]:
            if "diploid" in calldata[patient][sample]:
                dploidy = calldata[patient][sample]["diploid"][0]
                dmatch = "False"
                #This is where we'd add in something about 'the 2N peak was pretty wide' if we got that information.
                if abs(dploidy-2) < 0.3:
                    dmatch = "Two"
                else:
                    for flow in flowlist:
                        if abs(flow-dploidy) < 0.2:
                            dmatch = "True"
                            break
            tploidy = "None"
            if "eight" in calldata[patient][sample]:
                tploidy = calldata[patient][sample]["eight"][0]
            elif "tetraploid" in calldata[patient][sample]:
                tploidy = calldata[patient][sample]["tetraploid"][0]
            if tploidy != "None":
                tmatch = "False"
                for flow in flowlist:
                    if abs(flow-tploidy) < 0.2:
                        tmatch = "True"
                        break
    return (dmatch, tmatch)

def getBetterAccuracy(patient, sample, calldata):
    better_accuracy = "Unknown"
    if patient in calldata and sample in calldata[patient]:
        if "diploid" in calldata[patient][sample]:
            dacc = calldata[patient][sample]["diploid"][2]
            tacc = 0
            whichtet = "Tetraploid"
            if "tetraploid" in calldata[patient][sample]:
                tacc = calldata[patient][sample]["tetraploid"][2]
            if "eight" in calldata[patient]:
                whichtet = "Eight"
                tacc = calldata[patient][sample]["eight"][2]
            if tacc==0:
                better_accuracy = "Diploid only"
            else:
                if abs(dacc-tacc) < .03:
                    better_accuracy = "Neither"
                elif dacc > tacc:
                    better_accuracy = "Diploid"
                else:
                    better_accuracy = whichtet
        elif "eight" in calldata[patient][sample]:
            better_accuracy = "Eight only"
        elif "tetraploid" in calldata[patient][sample]:
            better_accuracy = "Tetraploid only"
    return better_accuracy

flowdata = readFlowSummary(flow_summary)
calldata = readCallSummary(call_summary)
outdata = open(outfile, "w")
for line in open(initcall, "r"):
    if "Patient" in line:
        writeNewLine(outdata, line, ("Flow ratio", "Close diploid flow?", "Close tetraploid flow?", "Better accuracy?"))
        continue
    line = line.replace(":", "::")
    (patient, sample, tomcall, vafcat) = line.split('\t')
    dmatch, tmatch = getDtMatches(patient, sample, flowdata, calldata)
    better_accuracy = getBetterAccuracy(patient, sample, calldata)
    
    writeNewLine(outdata, line, (flowdata[patient][0], dmatch, tmatch, better_accuracy))
    
outdata.close()

