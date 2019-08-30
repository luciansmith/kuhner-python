#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 11:25:43 2019

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy
from ete3 import Tree

import numpy
import math
import matplotlib.pyplot as plt
import csv

import lucianSNPLibrary as lsl

jabbafile = "JaBbA_SegmentedCopyNumber_Xiaotong/20190529_JaBbASegmentedAbsoluteCopyNumber_Xiaotong.txt"
CNdir = "noninteger_processed_CNs/"
eafile = "JaBbA_SegmentedCopyNumber_Xiaotong/EA_Drivers_Gene_List.txt"

CNfiles = []
for __, _, files in walk(CNdir):
    CNfiles += files
    
somepatientsonly = False
somepatients = ["74"]

def findHighCopyGenesFor(genelist, cncalls):
    retlist = {}
    for gene in genelist:
        retlist[gene] = []
        (chr, start, end) = genelist[gene]
        for patient in cncalls:
            for sample in cncalls[patient]:
#                if not sample=="24714":
#                    continue
                if chr in cncalls[patient][sample]:
                    totaloverlap = 0
                    for seg in cncalls[patient][sample][chr]:
                        if start > seg[1]:
                            continue
                        if end < seg[0]:
                            continue
                        overlapstart = max(start, seg[0])
                        overlapend = min(end, seg[1])
                        totaloverlap += overlapend-overlapstart
#                        print(gene)
#                        print("Overlap:", str(overlapstart), str(overlapend), "from", str(seg), str(start), str(end))
                        if totaloverlap > (end - start)/2:
                            retlist[gene].append((patient, sample))
                            continue
#                            print("Yes")
#                        else:
#                            print("No")
    return retlist

def writeFileFor(file, genelist):
    file.write("Gene\tHigh levels in:\n")
    for gene in genelist:
        file.write(gene)
        for (patient, sample) in genelist[gene]:
            file.write("\t" + patient + "_" + sample)
        file.write("\n")
    file.close()


jabbaCNs = {}
for line in open(jabbafile, "r"):
    if "track" in line:
        continue
    (id, chr, start, end, __, cn) = line.rstrip().split()
    (patient, sample) = id.split("-")[0:2]
    if somepatientsonly and patient not in somepatients:
        continue
    ploidy = lsl.getBestPloidyFor(patient, sample, challenge=False)
    cn = int(cn)
    
    if cn<7:
        continue
    start = int(start)
    end = int(end)
    if patient not in jabbaCNs:
        jabbaCNs[patient] = {}
    if sample not in jabbaCNs[patient]:
        jabbaCNs[patient][sample] = {}
    if chr not in jabbaCNs[patient][sample]:
        jabbaCNs[patient][sample][chr] = []
    jabbaCNs[patient][sample][chr].append([start, end, cn])

pascatCNs = {}
for file in CNfiles:
    (patient, sample, __, ploidy) = file.split("_")[0:4]
    if somepatientsonly and patient not in somepatients:
        continue
    bestploidy = lsl.getBestPloidyFor(patient, sample, challenge=False)
    if not bestploidy == ploidy:
        continue
    print("Getting sample", sample, "ploidy", ploidy)
    if patient not in pascatCNs:
        pascatCNs[patient] = {}
    for line in open(CNdir + file, "r"):
        if "patient" in line:
            continue
        (__, ___, chr, start, end, __, __, intA, intB) =line.rstrip().split()
        try:
            intA = int(intA)
            intB = int(intB)
        except:
            continue
        if ploidy=="tetraploid":
            if intA+intB < 7:
                continue
        elif intA+intB < 7:
            continue
        start = int(start)
        end = int(end)
        if sample not in pascatCNs[patient]:
            pascatCNs[patient][sample] = {}
        if chr not in pascatCNs[patient][sample]:
            pascatCNs[patient][sample][chr] = []
        pascatCNs[patient][sample][chr].append([start, end, intA+intB])

ea_drivers = {}
for line in open(eafile, "r"):
    if "Gene" in line:
        continue
    (gene, chr, start, end) = line.rstrip().split()
    try:
        int(chr)
    except:
        continue
    ea_drivers[gene] = (chr, int(start), int(end))

jlist = findHighCopyGenesFor(ea_drivers, jabbaCNs)
plist = findHighCopyGenesFor(ea_drivers, pascatCNs)

jfile = open("jabba_list.txt", "w")
writeFileFor(jfile, jlist)

pfile = open("pascat_list.txt", "w")
writeFileFor(pfile, plist)

        
        
