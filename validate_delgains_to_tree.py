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

tag = "g500_better_ploidy/"
delgain = "joint_delsandgains_"
delgaindir = delgain + tag

treedir = "phylip_TS_analysis/"

outdir = "delgains_validated/"

if not(path.isdir(outdir)):
    mkdir(outdir)


delgainlist = []
for (_, _, f) in walk(delgaindir):
    delgainlist += f

dgs = {}
for dgfilename in delgainlist:
    patient = dgfilename.split('_')[0]
    isA = "A"
    if "allB" in dgfilename:
        isA = "B"
    if patient not in dgs:
        dgs[patient] = {}
    dgfile = open(delgaindir + dgfilename, "r")
    samples = []
    for line in dgfile:
        lvec = line.rstrip().split()
        if "chr" in line:
            samples = lvec[4:]
            continue
        segid = lvec[0:4]
        segid.append(isA)
        segid = tuple(segid)
        dgs[patient][segid] = {}
        for i in range(len(samples)):
            sample = samples[i]
            dgs[patient][segid][sample] = lvec[i+4]
    dgfile.close()

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

allvalid = open(outdir + "all_valid.txt", "w")
allinvalid = open(outdir + "all_invalid.txt", "w")
allvalid.write  ("chr\tstart\tend\tcopynumber\tA or B\tValid\tnTrue\tnFalse")
allinvalid.write("chr\tstart\tend\tcopynumber\tA or B\tValid\tnTrue\tnFalse")
validlengths = {}
invalidlengths = {}
for i in range(16):
    if i==1:
        continue
    validlengths[i] = []
    invalidlengths[i] = []
for patient in dgs:
    treefilename = treedir + patient + "_outtree.txt"
    if not path.exists(treefilename):
        continue
    outfile = open(outdir + patient + "_delgain_validated.txt", "w")
    outfile.write("chr\tstart\tend\tcopynumber\tA or B\tValid\tnTrue\tnFalse")
    samples = list(dgs[patient][list(dgs[patient].keys())[0]].keys())
    for sample in samples:
        outfile.write("\t" + sample)
    outfile.write("\n")
    tree = ete3.Tree(treefilename)
    validate = {}
    setids = list(dgs[patient].keys())
    setids.sort()
    for segid in setids:
        truesamples = []
        falsesamples = []
        for sample in dgs[patient][segid]:
            if dgs[patient][segid][sample] == "True":
                truesamples.append(sample)
            elif dgs[patient][segid][sample] == "False":
                falsesamples.append(sample)
        valid_out = validateSamples(truesamples, falsesamples, tree)
        for seg in segid:
            outfile.write(seg + "\t")
        for validity in valid_out:
            outfile.write(str(validity) + "\t")
        for sample in samples:
            outfile.write(dgs[patient][segid][sample] + "\t")
        outfile.write("\n")
        if valid_out[0] == "Valid":
            for seg in segid:
                allvalid.write(seg + "\t")
            for validity in valid_out:
                allvalid.write(str(validity) + "\t")
            length = int(segid[2]) - int(segid[1])
            validlengths[int(segid[3])].append(length)
            allvalid.write(str(length) + "\n")
        if valid_out[0] == "Invalid":
            for seg in segid:
                allinvalid.write(seg + "\t")
            for validity in valid_out:
                allinvalid.write(str(validity) + "\t")
            length = int(segid[2]) - int(segid[1])
            invalidlengths[int(segid[3])].append(length)
            allinvalid.write(str(length) + "\n")
            
    outfile.close()

allvalid.close()
allinvalid.close()

for n in range(8):
    if n==1:
        continue
    if len(validlengths[n]) > 10:
        print("Valid lengths for copy number call", str(n))
        print("  number of calls:", str(len(validlengths[n])))
        print("  mean: ", str(np.mean(validlengths[n])))
        print("  median: ", str(np.median(validlengths[n])))
        print("  stdev: ", str(np.std(validlengths[n])))
        if len(validlengths[n]) > 30:
            x = lsl.createPrintAndSaveHistogram(validlengths[n], "", 10000)
    if len(invalidlengths[n]) > 10:
        print("Invalid lengths for copy number call", str(n))
        print("  number of calls:", str(len(invalidlengths[n])))
        print("  mean: ", str(np.mean(invalidlengths[n])))
        print("  median: ", str(np.median(invalidlengths[n])))
        print("  stdev: ", str(np.std(invalidlengths[n])))
        if len(invalidlengths[n]) > 30:
            y = lsl.createPrintAndSaveHistogram(invalidlengths[n], "", 10000)
