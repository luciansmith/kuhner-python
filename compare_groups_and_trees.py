# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:43:41 2018

@author: Lucian
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt
import csv
import ete3

import lucianSNPLibrary as lsl

highVAF = False
lowVAF = False

tag = ""
if highVAF:
    tag = "_highVAF"
if lowVAF:
    tag = "_lowVAF"

groupdir = "SNV_groups" + tag + "/"
treedir = "phylip_TS_analysis/"
outfilename = "patient_analysis.tsv"
allg_outfile = "all_groups.tsv"

CNVfile = "qc_consolidated.txt"
#CNVfile = "qc_xiaohong.txt"

#outdir = "SNV_CNV_tree_compare" + tag + "/"
#if not path.isdir(outdir):
#    mkdir(outdir)

groupfiles = []
for __, _, files in walk(groupdir):
    groupfiles += files

def callGroup(group, allsamples, tree):
    if len(group) == 1:
        return "Singleton"
    if len(group) == len(allsamples):
        return "Root"
    trueset = set()
    falseset = set()
    for branch in tree:
        if branch.name != "":
            for sample in group:
                if sample in branch.name:
                    trueset.add(branch)
            for sample in allsamples:
                if sample in group:
                    continue
                if sample in branch.name:
                    falseset.add(branch)
            if "blood" in branch.name:
                tree.set_outgroup(branch)
    trueroot = tree.get_common_ancestor(trueset)
    for fbranch in falseset:
        testset = trueset.copy()
        testset.add(fbranch)
        newroot = tree.get_common_ancestor(testset)
        if newroot == trueroot:
            return "Ungrouped"
    return "Grouped"


# Read in the SNV numbers
groupdata = {}
allsamples= {}
for gfile in groupfiles:
    if "all" in gfile or "patient" in gfile:
        continue
    patient = gfile.split("_")[0]
    groupdata[patient] = {}
    allsamples[patient] = set()
    for line in open(groupdir + gfile, "r"):
        if "Patient" in line:
            continue
        lvec = line.rstrip().split()
        assert(patient == lvec[0])
        count = int(lvec[1])
        perc = float(lvec[2])
        samples = tuple(lvec[3:])
        groupdata[patient][samples] = {}
        groupdata[patient][samples]["count"] = count
        groupdata[patient][samples]["percentage"] = perc
        groupdata[patient][samples]["cnv_count"] = 0
        groupdata[patient][samples]["cnv_percentage"] = 0.0
        for sample in samples:
            allsamples[patient].add(sample)

#Read in the CNV numbes
CNVs = {}
samplelists = {}
for line in open(CNVfile, "r"):
    if "chr" in line:
        continue
    lvec = line.rstrip().split()
    (patient, sample, __, __, __, nbaf, call, eval, accuracy, qc) = lvec
    if call == "1/1":
        continue
    segid = tuple(lvec[2:6])
    if patient not in CNVs:
        CNVs[patient] = {}
    if segid not in CNVs[patient]:
        CNVs[patient][segid] = {}
    if call not in CNVs[patient][segid]:
        CNVs[patient][segid][call] = []
    CNVs[patient][segid][call].append(sample)

#Count the CNVs by sample list
#nmulti = 0
#nmulti_singletons = 0
#nmulti_multis = 0
#nsingle = 0
CNVcounts = {}
for patient in CNVs:
    CNVtotal = 0
    for segid in CNVs[patient]:
        for call in CNVs[patient][segid]:
            samples = CNVs[patient][segid][call]
#        calls = list(CNVs[patient][segid].keys())
#        if len(calls) > 1:
#            nmulti += 1
#            nOne = 0
#            nTwoPlus = 0
#            for call in calls:
#                if len(CNVs[patient][segid][call])==1:
#                    nOne += 1
#                else:
#                    nTwoPlus += 1
#            if nTwoPlus == 0:
#                nmulti_singletons += 1
#            elif nOne == 0:
#                nmulti_multis += 1
#            
#            continue
#        samples = CNVs[patient][segid][calls[0]]
            samples.sort()
            samples = tuple(samples)
            if samples not in groupdata[patient]:
                groupdata[patient][samples] = {}
                groupdata[patient][samples]["count"] = 0
                groupdata[patient][samples]["percentage"] = 0.0
                groupdata[patient][samples]["cnv_count"] = 0
                groupdata[patient][samples]["cnv_percentage"] = 0.0
            groupdata[patient][samples]["cnv_count"] += 1
            CNVtotal += 1
#            nsingle += 1
    for samples in groupdata[patient]:
        groupdata[patient][samples]["cnv_percentage"] = groupdata[patient][samples]["cnv_count"]/CNVtotal

#print("Number of segments with a single call:", str(nsingle))
#print("Number of segments with multiple calls:", str(nmulti))
#print("Number of segments with multiple calls, all singletons:", str(nmulti_singletons))
#print("Number of segments with multiple calls, all multiples:", str(nmulti_multis))

#Now put the tree data in there, too:
for patient in groupdata:
    treefilename = treedir + patient + "_outtree.txt"
    if patient == "891":
        treefilename = treedir + patient + "a_outtree.txt"
    tree = ete3.Tree(treefilename)
    for samples in groupdata[patient]:
        groupdata[patient][samples]["matches_tree"] = callGroup(samples, allsamples[patient], tree)

#And finally, write out all of our information.
outfile = open(groupdir + allg_outfile, "w")
outfile.write("Patient\tMatches_tree\tCount\tPercentage\tCNV count\tCNV percentage\tSample1\tSample2\tSample3\tSample4\tSample5\tSample6\n")
for patient in groupdata:
    for samples in groupdata[patient]:
        outfile.write(patient)
        outfile.write("\t" + groupdata[patient][samples]["matches_tree"])
        outfile.write("\t" + str(groupdata[patient][samples]["count"]))
        outfile.write("\t" + str(groupdata[patient][samples]["percentage"]))
        outfile.write("\t" + str(groupdata[patient][samples]["cnv_count"]))
        outfile.write("\t" + str(groupdata[patient][samples]["cnv_percentage"]))
        for sample in samples:
            outfile.write("\t" + sample)
        outfile.write("\n")
outfile.close()

#Now do some analysis
has23GD = ["74", "279", "303", "391", "396", "450", "772", "997"]
types = ["Singleton", "Root", "Grouped", "Ungrouped"]
outfile = open(groupdir + outfilename, "w")
outfile.write("Patient\tnSNVmin\tnSNVmax\thas 2-3 GD")
for type in types:
    outfile.write("\t" + type + " counts")
    outfile.write("\t" + type + " total")
outfile.write("\tUngrouped potential subclone counts\tUngrouped potential scomubclone total\n")
for patient in groupdata:
    smallestSNVcount = 100000
    maxSNVcount = 0
    for samples in groupdata[patient]:
        SNVcount = groupdata[patient][samples]["count"]
        if groupdata[patient][samples]["matches_tree"] == "Grouped":
            if SNVcount < smallestSNVcount:
                smallestSNVcount = SNVcount
        if SNVcount > maxSNVcount:
            maxSNVcount = SNVcount
    possibleSubcloneThreshold = smallestSNVcount*3/4
    CNVcounts = {}
    for match in types:
        CNVcounts[match] = []
    CNVcounts["subclones"] = []
    for samples in groupdata[patient]:
        theseSamples = groupdata[patient][samples]
        CNVcount = theseSamples["cnv_count"]
        if CNVcount ==0:
            continue
        if theseSamples["matches_tree"] == "Ungrouped" and theseSamples["count"] >= possibleSubcloneThreshold:
            CNVcounts["subclones"].append(CNVcount)
        else:
            CNVcounts[theseSamples["matches_tree"]].append(CNVcount)
    outfile.write(patient)
    outfile.write("\t" + str(possibleSubcloneThreshold))
    outfile.write("\t" + str(maxSNVcount))
    outfile.write("\t" + str(patient in has23GD))
    for type in types:
        outfile.write("\t")
        for num in CNVcounts[type]:
            outfile.write("//" + str(num))
        outfile.write("\t" + str(sum(CNVcounts[type])))
    outfile.write("\t")
    for num in CNVcounts["subclones"]:
        outfile.write("//" + str(num))
    outfile.write("\t" + str(sum(CNVcounts["subclones"])))
    outfile.write("\n")
outfile.close()
    












