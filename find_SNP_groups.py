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

import lucianSNPLibrary as lsl

highVAF = False
lowVAF = False

tag = ""
if highVAF:
    tag = "_highVAF"
if lowVAF:
    tag = "_lowVAF"

mutation_file = "snv_plus_indels.twoPlus.20181030.csv"
outdir = "SNV_groups" + tag + "/"

if not path.isdir(outdir):
    mkdir(outdir)

def getPatientSampleMap():
    patientSampleMap = {}
    samplePatientMap = {}
    callfile = open("calling_evidence.tsv", "r")
    for line in callfile:
        if "Patient" in line:
            continue
        (patient, sample) = line.rstrip().split()[0:2]
        patientSampleMap[sample] = patient
        if patient not in samplePatientMap:
            samplePatientMap[patient] = []
        samplePatientMap[patient].append(sample)
    return patientSampleMap, samplePatientMap

patientSampleMap, samplePatientMap = getPatientSampleMap()

mutations = {}
with open(mutation_file, 'r') as csvfile:
    for lvec in csv.reader(csvfile):
        if "DNANum" in lvec[0]:
            continue
        (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
        refcnt = int(lvec[47])
        mutcnt = int(lvec[48])
        VAF = mutcnt/(refcnt+mutcnt)
        #print(VAF)
        if highVAF and VAF < 0.25:
            continue
        if lowVAF and VAF >= 0.25:
            continue
        if (is_snv=="f"):
            continue
        if (is_2p=="f"):
            continue
    #    if ("N" in sample):
    #        continue
        chr = int(chr)
        pos = int(pos)
        patient = patientSampleMap[sample]
        if patient not in mutations:
            mutations[patient] = {}
        if chr not in mutations[patient]:
            mutations[patient][chr] = {}
        if (pos, ref, alt) not in mutations[patient][chr]:
            mutations[patient][chr][(pos, ref, alt)] = set()
        mutations[patient][chr][(pos, ref, alt)].add(sample)

for patient in mutations:
    groupcounts = {}
    total = 0
    for chr in mutations[patient]:
        for posrefalt in mutations[patient][chr]:
            sampleset = mutations[patient][chr][posrefalt]
            sampleset = list(sampleset)
            sampleset.sort()
            sampleset = tuple(sampleset)
            if sampleset not in groupcounts:
                groupcounts[sampleset] = 0
            groupcounts[sampleset] += 1
            total += 1
    countfile = open(outdir + patient + "_SNV_counts.tsv", "w")
    countfile.write("Patient\tCount\tPercentage\tSamples\n")
    for sampleset in groupcounts:
        countfile.write(patient)
        countfile.write("\t" + str(groupcounts[sampleset]))
        countfile.write("\t" + str(groupcounts[sampleset]/total))
        for sample in sampleset:
            countfile.write("\t" + sample)
        countfile.write("\n")
    countfile.close()

