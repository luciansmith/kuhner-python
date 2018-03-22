#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 14:36:23 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os.path import isfile
from os import mkdir
from shutil import copy2 as copy
from copy import deepcopy

import numpy
import math
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import lucianSNPLibrary as lsl


indeldir = "/media/mkkuhner/Seagate Backup Plus Drive/wgs/all_wgs_SNV_Indel_txt_files/"
#indeldir = "all_wgs_SNV_Indel_txt_files/"
nonintcall_dir = "noninteger_processed_CNs/"
VAFpngdir = "VAF_pngs/"
VAF1v2dir = "VAF_1v2_hist/"

if not path.isdir(VAFpngdir):
    mkdir(VAFpngdir)

if not path.isdir(VAF1v2dir):
    mkdir(VAF1v2dir)

def getCallFiles():
    ret = {}
    files = []
    for (__, __, f) in walk(nonintcall_dir):
        files += f
    for f in files:
        if "nonint_CN" not in f:
            continue
        (patient, sample, gamma, ploidy, __, __) = f.split("_")
        ret[(patient, sample, ploidy)] = f
    return ret

def getCallsFor(patient, sample, files):
    calls = {}
    for ploidy in ("diploid", "tetraploid"):
        calls[ploidy] = {}
        if (patient, sample, ploidy) not in files:
            continue
        callfile = open(nonintcall_dir + files[(patient, sample, ploidy)], "r")
        for line in callfile:
            if "patient" in line:
                continue
            (__, __, chr, start, end, rawA, rawB, intA, intB) = line.split()
            if intA=="NA" or intB=="NA":
                continue
            if chr not in calls[ploidy]:
                calls[ploidy][chr] = []
            calls[ploidy][chr].append((int(start), int(end), int(intA) + int(intB)))
    return calls

def isOneVsTwo(calls, chr, pos):
    if "diploid" not in calls or "tetraploid" not in calls:
        return False
    if chr not in calls["diploid"] or chr not in calls["tetraploid"]:
        return False
    dipval = -1
    tetval = -1
    for seg in calls["diploid"][chr]:
        if pos >= seg[0] and pos <= seg[1]:
            dipval = seg[2]
            break
    if dipval != 1:
        return False
    for seg in calls["tetraploid"][chr]:
        if pos >= seg[0] and pos <= seg[1]:
            tetval = seg[2]
            break
    if tetval > 1:
        return True
    return False

callfiles = getCallFiles()
numpoints = open(VAF1v2dir + "summary.tsv", "w")
numpoints.write("Patient")
numpoints.write("\tSample")
numpoints.write("\tAll VAFS")
numpoints.write("\t1v2 VAFs")
numpoints.write("\n")

files = []
for (__, __, f) in walk(indeldir):
    files += f
for f in files:
    if "annotated" not in f:
        continue
    (patient, sample) = f.split("-")[0:2]
    calls = getCallsFor(patient, sample, callfiles)
    allVAFs = []
    onevtwo_VAFs = []
    indelfile = open(indeldir + f, "r")
    for line in indelfile:
        if "CHROM" in line:
            continue
        lvec = line.split()
        chrom = lvec[0]
        if chrom=="X" or chrom=="Y":
            continue
        pos = int(lvec[1])
        called = lvec[5]
        talt = int(lvec[30])
        tref = int(lvec[31])
        if talt < 0 or tref < 0:
            continue
        VAF = talt / (talt+tref)
        allVAFs.append(VAF)
        if(isOneVsTwo(calls, chrom, int(pos))):
            onevtwo_VAFs.append(VAF)
    #print("All VAFs for patient", patient, "sample", sample, ":")
    #hist = lsl.createPrintAndSaveHistogram(allVAFs, VAFpngdir + patient + "_" + sample + "_VAF_hist", 0.001, xdata="VAF", show=False, savefig=True)
    print("VAFs for positions called 01 in diploid but more in tetraploid,", patient, "sample", sample, ":")
    hist = lsl.createPrintAndSaveHistogram(onevtwo_VAFs, VAF1v2dir + patient + "_" + sample + "_1v2_VAF_hist", 0.001, xdata="VAF", show=False, savefig=True, axis=(0, 1.1, 0))
    numpoints.write(patient)
    numpoints.write("\t" + sample)
    numpoints.write("\t" + str(len(allVAFs)))
    numpoints.write("\t" + str(len(onevtwo_VAFs)))
    numpoints.write("\n")

numpoints.close()
