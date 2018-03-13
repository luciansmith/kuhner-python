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


#indeldir = "/media/mkkuhner/Seagate Backup Plus Drive/wgs/all_wgs_SNV_Indel_txt_files/"
indeldir = "all_wgs_SNV_Indel_txt_files/"
VAFpngdir = "VAF_histograms/"


if not path.isdir(VAFpngdir):
    mkdir(VAFpngdir)

files = []
for (__, __, f) in walk(indeldir):
    files += f
for f in files:
    if "annotated" not in f:
        continue
    (patient, sample) = f.split("-")[0:2]
    allVAFs = []
    twoplus_VAFs = []
    indelfile = open(indeldir + f, "r")
    for line in indelfile:
        if "CHROM" in line:
            continue
        lvec = line.split()
        chrom = lvec[0]
        if chrom=="X" or chrom=="Y":
            continue
        called = lvec[5]
        talt = int(lvec[30])
        tref = int(lvec[31])
        if talt < 0 or tref < 0:
            continue
        VAF = talt / (talt+tref)
        allVAFs.append(VAF)
        if len(called.split("-")) > 1:
            twoplus_VAFs.append(VAF)
    print("All VAFs for patient", patient, "sample", sample, ":")
    hist = lsl.createPrintAndSaveHistogram(allVAFs, VAFpngdir + patient + "_" + sample + "_VAF_hist", 0.001, xdata="VAF", show=False, savefig=True)
    peaks = lsl.getPeaks(hist)
#    print("Only VAFs from 2+ callers for patient", patient, "sample", sample, ":")
#    lsl.createPrintAndSaveHistogram(twoplus_VAFs, "", 0.001, xdata="VAF")
