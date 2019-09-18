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

#import lucianSNPLibrary as lsl

groupdir = "SNV_groups/"
treedir = "phylip_TS_analysis/"
SVfile = "SV_events.txt"
patientfile = "patient_analysis_SVs.tsv"
allg_outfile = "all_groups_SVs.tsv"

#outdir = "SNV_SV_tree_compare" + tag + "/"
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
    if "all" in gfile or "patient" in gfile or "SV_" in gfile:
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
        groupdata[patient][samples]["SV_count"] = 0
        for sample in samples:
            allsamples[patient].add(sample)

#Read in the SV numbes
SVs = {}
samplelists = {}
for line in open(SVfile, "r"):
    if "chr" in line:
        continue
    lvec = line.rstrip().split()
    (__, __, patient, sample, type, ch1, start1, end1, __, ch2, start2, end2, __, __, __, __) = lvec
    svid = (type, ch1, start1, end1, ch2, start2, end2)
    if patient not in SVs:
        SVs[patient] = {}
    if svid not in SVs[patient]:
        SVs[patient][svid] = set()
    SVs[patient][svid].add(sample)

#Count the SVs by sample list
nmulti = 0
nmulti_singletons = 0
nmulti_multis = 0
nsingle = 0
SVcounts = {}
for patient in SVs:
    SVtotal = 0
    for segid in SVs[patient]:
        samples = list(SVs[patient][segid])
        samples.sort()
        samples = tuple(samples)
        if samples not in groupdata[patient]:
            groupdata[patient][samples] = {}
            groupdata[patient][samples]["count"] = 0
            groupdata[patient][samples]["percentage"] = 0.0
            groupdata[patient][samples]["SV_count"] = 0
        groupdata[patient][samples]["SV_count"] += 1
        SVtotal += 1
        nsingle += 1
    for samples in groupdata[patient]:
        groupdata[patient][samples]["SV_percentage"] = groupdata[patient][samples]["SV_count"]/SVtotal

print("Number of segments with a single call:", str(nsingle))
print("Number of segments with multiple calls:", str(nmulti))
print("Number of segments with multiple calls, all singletons:", str(nmulti_singletons))
print("Number of segments with multiple calls, all multiples:", str(nmulti_multis))

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
outfile.write("Patient\tMatches_tree\tCount\tPercentage\tSV count\tSV percentage\tSample1\tSample2\tSample3\tSample4\tSample5\tSample6\n")
for patient in groupdata:
    for samples in groupdata[patient]:
        outfile.write(patient)
        outfile.write("\t" + groupdata[patient][samples]["matches_tree"])
        outfile.write("\t" + str(groupdata[patient][samples]["count"]))
        outfile.write("\t" + str(groupdata[patient][samples]["percentage"]))
        outfile.write("\t" + str(groupdata[patient][samples]["SV_count"]))
        outfile.write("\t" + str(groupdata[patient][samples]["SV_percentage"]))
        for sample in samples:
            outfile.write("\t" + sample)
        outfile.write("\n")
outfile.close()

#Now do some analysis
has23GD = ["74", "279", "303", "391", "396", "450", "772", "997"]
types = ["Singleton", "Root", "Grouped", "Ungrouped"]
outfile = open(groupdir + patientfile, "w")
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
    SVcounts = {}
    for match in types:
        SVcounts[match] = []
    SVcounts["subclones"] = []
    for samples in groupdata[patient]:
        theseSamples = groupdata[patient][samples]
        SVcount = theseSamples["SV_count"]
        if SVcount ==0:
            continue
        if theseSamples["matches_tree"] == "Ungrouped" and theseSamples["count"] >= possibleSubcloneThreshold:
            SVcounts["subclones"].append(SVcount)
        else:
            SVcounts[theseSamples["matches_tree"]].append(SVcount)
    outfile.write(patient)
    outfile.write("\t" + str(possibleSubcloneThreshold))
    outfile.write("\t" + str(maxSNVcount))
    outfile.write("\t" + str(patient in has23GD))
    for type in types:
        outfile.write("\t")
        for num in SVcounts[type]:
            outfile.write("//" + str(num))
        outfile.write("\t" + str(sum(SVcounts[type])))
    outfile.write("\t")
    for num in SVcounts["subclones"]:
        outfile.write("//" + str(num))
    outfile.write("\t" + str(sum(SVcounts["subclones"])))
    outfile.write("\n")
outfile.close()
    












