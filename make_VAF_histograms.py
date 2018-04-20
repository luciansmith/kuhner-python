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


runLocally = False
justdip = False
just_intersection = False
onlysomepatients = False
somepatients = ["954"]

import matplotlib

if runLocally:
    indeldir = "all_wgs_SNV_Indel_txt_files/"
else:
    matplotlib.use('Agg')
    indeldir = "/media/mkkuhner/Seagate Backup Plus Drive/wgs/all_wgs_SNV_Indel_txt_files/"

import matplotlib.pyplot as plt
import lucianSNPLibrary as lsl

nonintcall_dir = "noninteger_processed_CNs/"
VAFpngdir = "VAF_pngs/"
VAF1v2dir = "VAF_1v2_hist/"
VAF2v4dir = "VAF_2v4_hist/"
XVAFdir = "XVAF_hist/"
dipCNLOHdir = "dip_CNLOH_hist/"
intersection_dir = "2v4_intersection/"

if not path.isdir(VAFpngdir):
    mkdir(VAFpngdir)

if not path.isdir(VAF1v2dir):
    mkdir(VAF1v2dir)

if not path.isdir(VAF2v4dir):
    mkdir(VAF2v4dir)

if not path.isdir(XVAFdir):
    mkdir(XVAFdir)

if not path.isdir(dipCNLOHdir):
    mkdir(dipCNLOHdir)

if not path.isdir(intersection_dir):
    mkdir(intersection_dir)

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

def getDuplicates(files):
    duplicates = []
    sumdata = {}
    sumfile = open("xiaocompare_summary.tsv", "r")
    for line in sumfile:
        if "Patient" in line:
            continue
        (patient, sample, gamma, ploidy, ploidyval) = line.split()[0:5]
        if (patient, sample) not in sumdata:
            sumdata[patient, sample] = {}
        sumdata[patient, sample][ploidy] = ploidyval
    for pat_sam in sumdata:
        if "diploid" in sumdata[pat_sam] and "tetraploid" in sumdata[pat_sam]:
            if sumdata[pat_sam]["diploid"] == sumdata[pat_sam]["tetraploid"]:
                duplicates.append(pat_sam)
    return duplicates

def getCallsFor(patient, sample, files, duplicates):
    calls = {}
    for ploidy in ("diploid", "tetraploid"):
        calls[ploidy] = {}
        if (patient, sample, ploidy) not in files:
            continue
        callfile = open(nonintcall_dir + files[(patient, sample, ploidy)], "r")
        if ploidy=="tetraploid":
            if (patient, sample, "eight") in files:
                callfile = open(nonintcall_dir + files[(patient, sample, "eight")], "r")
            elif (patient, sample) in duplicates:
                continue
        for line in callfile:
            if "patient" in line:
                continue
            (__, __, chr, start, end, rawA, rawB, intA, intB) = line.split()
            if intA=="NA" or intB=="NA":
                continue
            if chr not in calls[ploidy]:
                calls[ploidy][chr] = []
            calls[ploidy][chr].append((int(start), int(end), int(intA), int(intB)))
    return calls

def getTypeComparison(calls, chr, pos):
    dip_isloss = False
    dip_isdip = False
    tet_istwoplus = False
    tet_istet = False
    if "diploid" not in calls or chr not in calls["diploid"]:
        dip_isloss = True
        dip_isdip = True
    else:
        for seg in calls["diploid"][chr]:
            if pos >= seg[0] and pos <= seg[1]:
                dip_isloss = seg[2] + seg[3] == 1
                dip_isdip = seg[2]==1 and seg[3]==1
                if (seg[2]==0 and seg[3]==2) or (seg[3]==0 and seg[2]==2):
                    return "dip_CNLOH"
                break
        if not dip_isloss and not dip_isdip:
            return "None"
    if "tetraploid" not in calls or chr not in calls["tetraploid"]:
        tet_istwoplus = True
        tet_istet = True
    else:
        for seg in calls["tetraploid"][chr]:
            if pos >= seg[0] and pos <= seg[1]:
                tet_istwoplus = seg[2] + seg[3] >= 2 and seg[2] != seg[3]
                tet_istet = seg[2]==2 and seg[3]==2
                break
    if dip_isdip and tet_istet:
        return "2v4"
    if dip_isloss and tet_istwoplus:
        return "1v2"
    return "None"

callfiles = getCallFiles()
duplicates = getDuplicates(callfiles)
if not justdip:
    numpoints = open(VAF1v2dir + "summary.tsv", "w")
    numpoints.write("Patient")
    numpoints.write("\tSample")
    numpoints.write("\tAll VAFS")
    numpoints.write("\t1v2 VAFs")
    numpoints.write("\t1v2 VAFs > 0.55")
    numpoints.write("\t1v2 VAFs > 0.8")
    numpoints.write("\tX VAFs")
    numpoints.write("\tX VAFs > 0.55")
    numpoints.write("\tX VAFs > 0.8")
    numpoints.write("\t2v4 VAFs")
    numpoints.write("\t2v4 VAFs between 0.22 and 0.28")
    numpoints.write("\n")

CNLOH_test = open(VAF1v2dir + "dip_CNLOH.tsv", "w")
CNLOH_test.write("Patient")
CNLOH_test.write("\tSample")
CNLOH_test.write("\tAll VAFS")
CNLOH_test.write("\tdiploid CNLOH VAFs")
CNLOH_test.write("\tdiploid CNLOH VAFs > 0.55")
CNLOH_test.write("\tdiploid CNLOH VAFs > 0.8")
CNLOH_test.write("\n")

files = []
for (__, __, f) in walk(indeldir):
    files += f
for f in files:
    if "annotated" not in f:
        continue
    (patient, sample) = f.split("-")[0:2]
    if onlysomepatients and patient not in somepatients:
        continue
    calls = getCallsFor(patient, sample, callfiles, duplicates)
    allVAFs = []
    onevtwo_VAFs = []
    twovfour_VAFs = []
    dip_CNLOH_VAFs = []
    X_VAFs = []
    indelfile = open(indeldir + f, "r")
    print("Reading", f)
    intersection_file = open(intersection_dir + patient + "_" + sample + "_2v4_positions.tsv", "w")
    intersection_file.write("Patient")
    intersection_file.write("\tSample")
    intersection_file.write("\tChr")
    intersection_file.write("\tPosition")
    intersection_file.write("\n")
    for line in indelfile:
        if "CHROM" in line:
            continue
        lvec = line.split()
        chrom = lvec[0]
        pos = int(lvec[1])
        called = lvec[5]
        talt = int(lvec[30])
        tref = int(lvec[31])
        if talt < 0 or tref < 0:
            continue
        VAF = talt / (talt+tref)
        if chrom=="X":
            X_VAFs.append(VAF)
            continue
        if chrom=="Y":
            continue
        allVAFs.append(VAF)
        compare = getTypeComparison(calls, chrom, int(pos))
        if compare=="1v2":
            onevtwo_VAFs.append(VAF)
        elif compare=="2v4":
            twovfour_VAFs.append(VAF)
            intersection_file.write(patient)
            intersection_file.write("\t" + sample)
            intersection_file.write("\t" + chrom)
            intersection_file.write("\t" + str(pos))
            intersection_file.write("\n")
        elif compare=="dip_CNLOH":
            dip_CNLOH_VAFs.append(VAF)
    intersection_file.close()
    if just_intersection:
        continue
    if not justdip:
        #print("All VAFs for patient", patient, "sample", sample, ":")
        hist = lsl.createPrintAndSaveHistogram(allVAFs, VAFpngdir + patient + "_" + sample + "_VAF_hist", 0.001, xdata="VAF", show=runLocally, savefig=True, axis=(0, 1.1, 0))
        #print("VAFs for positions called 1/1 in diploid and 2/2 in tetraploid,", patient, "sample", sample, ":")
        hist = lsl.createPrintAndSaveHistogram(twovfour_VAFs, VAF2v4dir + patient + "_" + sample + "_2v4_VAF_hist", 0.001, xdata="VAF", show=runLocally, savefig=True, axis=(0, 1.1, 0))
        #print("VAFs for positions called 01 in diploid but more in tetraploid,", patient, "sample", sample, ":")
        hist = lsl.createPrintAndSaveHistogram(onevtwo_VAFs, VAF1v2dir + patient + "_" + sample + "_1v2_VAF_hist", 0.001, xdata="VAF", show=runLocally, savefig=True, axis=(0, 1.1, 0))
        #print("VAFs for positions on the X chromosome,", patient, "sample", sample, ":")
        hist = lsl.createPrintAndSaveHistogram(X_VAFs, XVAFdir + patient + "_" + sample + "_X_VAF_hist", 0.001, xdata="VAF", show=runLocally, savefig=True, axis=(0, 1.1, 0))
        numpoints.write(patient)
        numpoints.write("\t" + sample)
        numpoints.write("\t" + str(len(allVAFs)))
        numpoints.write("\t" + str(len(onevtwo_VAFs)))
        numpoints.write("\t" + str(sum(1 for x in onevtwo_VAFs if x>0.55)))
        numpoints.write("\t" + str(sum(1 for x in onevtwo_VAFs if x>0.8)))
        numpoints.write("\t" + str(len(X_VAFs)))
        numpoints.write("\t" + str(sum(1 for x in X_VAFs if x>0.55)))
        numpoints.write("\t" + str(sum(1 for x in X_VAFs if x>0.8)))
        numpoints.write("\t" + str(len(twovfour_VAFs)))
        numpoints.write("\t" + str(sum(1 for x in twovfour_VAFs if x>0.22 and x<0.28)))
        
        numpoints.write("\n")

#    print("VAFs for CNLOH positions in the diploid version,", patient, "sample", sample, ":")
#    hist = lsl.createPrintAndSaveHistogram(dip_CNLOH_VAFs, dipCNLOHdir + patient + "_" + sample + "_dip_CNLOH_hist", 0.001, xdata="VAF", show=runLocally, savefig=False, axis=(0, 1.1, 0))
    CNLOH_test.write(patient)
    CNLOH_test.write("\t" + sample)
    CNLOH_test.write("\t" + str(len(allVAFs)))
    CNLOH_test.write("\t" + str(len(dip_CNLOH_VAFs)))
    CNLOH_test.write("\t" + str(sum(1 for x in dip_CNLOH_VAFs if x>0.55)))
    CNLOH_test.write("\t" + str(sum(1 for x in dip_CNLOH_VAFs if x>0.8)))
    CNLOH_test.write("\n")

if not justdip:
    numpoints.close()
CNLOH_test.close()
