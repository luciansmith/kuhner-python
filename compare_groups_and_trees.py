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
CNVfile = "qc_events.txt"

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
    if "all" in gfile:
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
nmulti = 0
nsingle = 0
CNVcounts = {}
for patient in CNVs:
    CNVtotal = 0
    for segid in CNVs[patient]:
        calls = list(CNVs[patient][segid].keys())
        if len(calls) > 1:
            nmulti += 1
            continue
        samples = CNVs[patient][segid][calls[0]]
        samples.sort()
        samples = tuple(samples)
        if samples not in groupdata[patient]:
            groupdata[patient][samples] = {}
            groupdata[patient][samples]["count"] = 0
            groupdata[patient][samples]["percentage"] = 0.0
            groupdata[patient][samples]["cnv_count"] = 0
        groupdata[patient][samples]["cnv_count"] += 1
        CNVtotal += 1
        nsingle += 1
    for samples in groupdata[patient]:
        groupdata[patient][samples]["cnv_percentage"] = groupdata[patient][samples]["cnv_count"]/CNVtotal

print("Number of segments with multiple calls:", str(nmulti))
print("Number of segments with a single call:", str(nsingle))

#Now put the tree data in there, too:
for patient in groupdata:
    treefilename = treedir + patient + "_outtree.txt"
    if patient == "891":
        treefilename = treedir + patient + "a_outtree.txt"
    tree = ete3.Tree(treefilename)
    for samples in groupdata[patient]:
        groupdata[patient][samples]["matches_tree"] = callGroup(samples, allsamples[patient], tree)

#And finally, write out all of our information.
outfile = open(groupdir + "all_groups.tsv", "w")
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