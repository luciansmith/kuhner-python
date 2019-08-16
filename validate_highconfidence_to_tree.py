#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:53:45 2019

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
import ete3
import numpy as np
import lucianSNPLibrary as lsl

allconf = True

maryfile = "qc_events.txt"

treedir = "phylip_TS_analysis/"

outdir = "highconf_validated/"
if (allconf):
    outdir = "allconf_validated/"

if not(path.isdir(outdir)):
    mkdir(outdir)


def validateSamples(truesamples, falsesamples, tree):
    trueset = set()
    falseset = set()
    for branch in tree:
        if branch.name != "":
            for sample in truesamples:
                if sample in branch.name:
                    trueset.add(branch)
            for sample in falsesamples:
                if sample in branch.name:
                    falseset.add(branch)
            if "blood" in branch.name:
                tree.set_outgroup(branch)
    truelen = len(trueset)
    falselen = len(falseset)
    if truelen == 0 or truelen==1 or falselen==0:
        return ["Unvalidatable", truelen, falselen]
    trueroot = tree.get_common_ancestor(trueset)
    for fbranch in falseset:
        assert(len(trueset)==truelen)
        testset = trueset.copy()
        testset.add(fbranch)
        newroot = tree.get_common_ancestor(testset)
        if newroot == trueroot:
            return ["Invalid", truelen, falselen]
    return ["Valid", truelen, falselen]
    #Now get the set of children of that root

def getBestCall(calllist):
    bestcall = ""
    nbestcalls = 0
    nbestfails = 100
    for call in calllist:
        nfails = 0
        for (sample, eval, accuracy, qc) in calllist[call]:
            if qc=="FAIL":
                nfails += 1
        if len(calllist[call]) > nbestcalls:
            bestcall = call
            nbestfails = nfails
            continue
        if len(calllist[call]) == nbestcalls and not(allconf):
            if nfails < nbestfails:
                bestcall = call
                nbestfails = nfails
                continue
    return bestcall

segs = {}
samplelists = {}
for line in open(maryfile, "r"):
    if "chr" in line:
        continue
    lvec = line.rstrip().split()
    (patient, sample, __, __, __, nbaf, call, eval, accuracy, qc) = lvec
    if call == "1/1":
        continue
    segid = tuple(lvec[2:6])
    if patient not in segs:
        segs[patient] = {}
        samplelists[patient] = set()
    if segid not in segs[patient]:
        segs[patient][segid] = {}
    if call not in segs[patient][segid]:
        segs[patient][segid][call] = []
    samplelists[patient].add(sample)
    segs[patient][segid][call].append([sample, eval, accuracy, qc])

allvalid = open(outdir + "all_valid.txt", "w")
allinvalid = open(outdir + "all_invalid.txt", "w")
allvalid.write  ("chr\tstart\tend\tcopynumber\tcall\tValid\tnTrue\tnFalse\n")
allinvalid.write("chr\tstart\tend\tcopynumber\tcall\tValid\tnTrue\tnFalse\n")
validlengths = {}
invalidlengths = {}
nBafValid = []
nBafInvalid = []

mostCalls = ["0/0", "0/1", "0/2", "1/2", "1/3", "2/2", "2/3"]
byPatientValidity = open(outdir + "byPatientValidity.tsv", "w")
byPatientValidity.write("Patient")
for call in mostCalls:
    byPatientValidity.write("\t" + call + "_valid")
    byPatientValidity.write("\t" + call + "_invalid")
byPatientValidity.write("\n")

for patient in segs:
    nCallsValid = {}
    nCallsInvalid = {}
    samples = samplelists[patient]
    treefilename = treedir + patient + "_outtree.txt"
    if not path.exists(treefilename):
        continue
    outfile = open(outdir + patient + "_highconf.txt", "w")
    outfile.write("chr\tstart\tend\tcopynumber\tcall\tValid\tnTrue\tnFalse")
#    for sample in samplelists[patient]:
#        outfile.write("\t" + sample)
    outfile.write("\n")
    tree = ete3.Tree(treefilename)
    validate = {}
    setids = list(segs[patient].keys())
    setids.sort()
    for segid in segs[patient]:
        bestcall = getBestCall(segs[patient][segid])
        truesamples = []
        falsesamples = []
        calledsamples = set()
        for call in segs[patient][segid]:
            for (sample, eval, accuracy, qc) in segs[patient][segid][call]:
                if call == bestcall:
                    if allconf or qc=="PASS":
                        truesamples.append(sample)
                calledsamples.add(sample)
        for sample in samples:
            if sample not in calledsamples:
                falsesamples.append(sample)
        valid_out = validateSamples(truesamples, falsesamples, tree)
        for seg in segid:
            outfile.write(seg + "\t")
        for validity in valid_out:
            outfile.write(str(validity) + "\t")
#        for sample in samplelists[patient]:
#            outfile.write(segs[patient][segid][sample] + "\t")
        outfile.write("\n")
        if valid_out[0] == "Valid":
            for seg in segid:
                allvalid.write(seg + "\t")
            for validity in valid_out:
                allvalid.write(str(validity) + "\t")
            length = int(segid[2]) - int(segid[1])
            if call not in validlengths:
                validlengths[bestcall] = []
            validlengths[bestcall].append(length)
            allvalid.write(str(length) + "\n")
            if call not in nCallsValid:
                nCallsValid[bestcall] = 0
            nCallsValid[bestcall] += 1
            nBafValid.append(int(segid[3]))
        if valid_out[0] == "Invalid":
            for seg in segid:
                allinvalid.write(seg + "\t")
            for validity in valid_out:
                allinvalid.write(str(validity) + "\t")
            length = int(segid[2]) - int(segid[1])
            if call not in invalidlengths:
                invalidlengths[bestcall] = []
            invalidlengths[bestcall].append(length)
            allinvalid.write(str(length) + "\n")
            if call not in nCallsInvalid:
                nCallsInvalid[bestcall] = 0
            nCallsInvalid[bestcall] += 1
            nBafInvalid.append(int(segid[3]))
    byPatientValidity.write(patient)
    for call in mostCalls:
        if call in nCallsValid:
            byPatientValidity.write("\t" + str(nCallsValid[call]))
        else:
            byPatientValidity.write("\t0")
        if call in nCallsInvalid:
            byPatientValidity.write("\t" + str(nCallsInvalid[call]))
        else:
            byPatientValidity.write("\t0")
    byPatientValidity.write("\n")

byPatientValidity.close()

allvalid.close()
allinvalid.close()

#for call in validlengths:
#    print("Valid lengths for copy number call", call)
#    print("  number of calls:", str(len(validlengths[call])))
#    print("  mean: ", str(np.mean(validlengths[call])))
#    print("  median: ", str(np.median(validlengths[call])))
#    print("  stdev: ", str(np.std(validlengths[call])))
#    if len(validlengths[call]) > 30:
#        x = lsl.createPrintAndSaveHistogram(validlengths[call], "", 10000)
#for call in invalidlengths:
#    print("Invalid lengths for copy number call", call)
#    print("  number of calls:", str(len(invalidlengths[call])))
#    print("  mean: ", str(np.mean(invalidlengths[call])))
#    print("  median: ", str(np.median(invalidlengths[call])))
#    print("  stdev: ", str(np.std(invalidlengths[call])))
#    if len(invalidlengths[call]) > 30:
#        y = lsl.createPrintAndSaveHistogram(invalidlengths[call], "", 10000)

for call in validlengths:
        print("Comparison for call", call, ":")
        print("Valid: ", str(len(validlengths[call])))
        if call in invalidlengths:
            print("Invalid: ", str(len(invalidlengths[call])))
            print("Percent valid:", str(len(validlengths[call])/(len(validlengths[call])+len(invalidlengths[call]))))
        else:
            print("No invalids at all")