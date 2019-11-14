# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:43:41 2018

@author: Lucian
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

highVAF = False
lowVAF = False
highVAFcutoff = 0.25
lowVAFcutoff = 0.25
nMutationsCutoff = 0

tag = ""
if highVAF:
    tag = "_highVAF"
if lowVAF:
    tag = "_lowVAF"


VAFdir = "VAFclusters_kanika/"
outdir = "SNV_groups_kanika" + tag + "/"

if not path.isdir(outdir):
    mkdir(outdir)

VAFfiles = []
for __, _, files in walk(VAFdir):
    VAFfiles += files

def makeTuple(in_set):
    ret = list(in_set)
    ret.sort()
    return tuple(ret)

allclusters = {}
nclusters = {}
for file in VAFfiles:
    patient, sample = file.split("_")[0:2]
    if patient not in allclusters:
        allclusters[patient] = []
        nclusters[patient] = {}
    clusters = []
    for line in open(VAFdir + file, "r"):
        lvec = line.rstrip().split("\t")
        if "Patient" in line:
            headers = lvec[5:]
            for cluster in headers:
                cluster = eval(cluster)
                cluster = set(cluster)
                clusters.append(cluster)
            continue
        vafs = lvec[5:]
        for n, vaf in enumerate(vafs):
            if vaf=="":
                continue
            vaf = float(vaf)
            if highVAF and vaf<=highVAFcutoff:
                continue
            if lowVAF and vaf > lowVAFcutoff:
                continue
            cluster = clusters[n]
            if cluster in allclusters[patient]:
                continue
            if makeTuple(cluster) not in nclusters[patient]:
                nclusters[patient][makeTuple(cluster)] = 0
            nclusters[patient][makeTuple(cluster)] += 1
    for cluster in clusters:
        allclusters[patient].append(cluster)


for patient in nclusters:
    groupcounts = {}
    total = 0
    for cluster in nclusters[patient]:
        total += nclusters[patient][makeTuple(cluster)]
    countfile = open(outdir + patient + "_SNV_counts.tsv", "w")
    countfile.write("Patient\tCount\tPercentage\tSamples\n")
    for cluster in nclusters[patient]:
        if nclusters[patient][makeTuple(cluster)] < nMutationsCutoff:
            continue
        countfile.write(patient)
        countfile.write("\t" + str(nclusters[patient][makeTuple(cluster)]))
        countfile.write("\t" + str(nclusters[patient][makeTuple(cluster)]/total))
        for sample in cluster:
            if "-" not in sample:
                countfile.write("\t" + sample)
        for sample in cluster:
            if "-" in sample:
                countfile.write("\t" + sample)
        countfile.write("\n")
    countfile.close()

