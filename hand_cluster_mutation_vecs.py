#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 15:45:51 2018

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

onlysomepatients = True
somepatients = ["160"]

VAFdir = "VAFcompare/"

def getPatientCNVs(patient):
    if patient=="160":
        return [(1,92472295, 112983230), (3, 44000109, 56994234), (5, 0, 190000000), (9, 0, 32764933), (17, 0, 22190573), (23, 0, 154000000)]
    return [(23, 0, 154000000)]

flist = []
for (_, _, f) in walk(VAFdir):
    flist += f
    
for f in flist:
    if "allsamples" not in f:
        continue
    patient = f.split("_")[0]
    patientCNVs = getPatientCNVs(patient)
    if onlysomepatients and patient not in somepatients:
        continue
    clusters = {}
    VAFfile = open(VAFdir + f, "r")
    samplevec = []
    nsamples = 0
    for line in VAFfile:
        linevec = line.rstrip().split("\t")
        nsamples = int((len(linevec)-4)/2)
        if "Patient" in line:
            samplevec = linevec[4:4+nsamples]
            continue
        (patient, chr, pos, alt) = linevec[0:4]
        chr = int(chr)
        pos = int(pos)
        VAFs = linevec[4:4+nsamples]
        presenceVec = tuple(linevec[4+nsamples:])
        if presenceVec not in clusters:
            clusters[presenceVec] = {}
            for sample in samplevec:
                clusters[presenceVec][sample] = {}
                clusters[presenceVec][sample]["alldip"] = []
                for CNV in patientCNVs:
                    clusters[presenceVec][sample][CNV] = []
        whichCNV = "alldip"
        for (cnv_chr, start, end) in patientCNVs:
            if not(chr==cnv_chr):
                continue
            if pos >= start and pos <= end:
                whichCNV = (cnv_chr, start, end)
                break
        for n in range(0,nsamples):
            sample = samplevec[n]
            VAF = float(VAFs[n])
            clusters[presenceVec][sample][whichCNV].append(VAF)
    
    
    print("Looking at clusters for patient", patient, "\nSamples are:", str(samplevec))
    for presenceVec in clusters:
        print("Looking at clusters for vector", str(presenceVec))
        for n in range(0,nsamples):
            if presenceVec[n] == "True":
                sample = samplevec[n]
                if len(clusters[presenceVec][sample]["alldip"]) < 100:
                    print("Too few mutations in this cluster for sample", sample)
                    continue
                print("Mutations in this vector for sample", sample)
                lsl.createPrintAndSaveHistogram(clusters[presenceVec][sample]["alldip"], "", 0.001)
    for CNV in patientCNVs:
        print("Looking at clusters for patient", patient, ", segment", str(CNV), "\nSamples are:", str(samplevec))
        for presenceVec in clusters:
            print("Looking at clusters for vector", str(presenceVec))
            for n in range(0,nsamples):
                if presenceVec[n] == "True":
                    sample = samplevec[n]
                    if len(clusters[presenceVec][sample][CNV]) < 10:
                        print("Too few mutations in this cluster for sample", sample)
                        continue
                    print("Mutations in this vector for sample", sample)
                    lsl.createPrintAndSaveHistogram(clusters[presenceVec][sample][CNV], "", 0.001)
        
        
