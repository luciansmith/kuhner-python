#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 12:49:23 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
import sys

import lucianSNPLibrary as lsl

use_baf = True
use_lengths = False
use_disjoin = False

if (use_baf and use_lengths) or (use_baf and use_disjoin) or (use_lengths and use_disjoin):
    print "You have to pick whether to run this over BAF values, length values, or disjoin values.\n"
    sys.exit()

use_min = False
nsamples_min = 0 #Arbitrary value: minimum number of samples we require
use_max = False
nsamples_max = 10000000 #Arbitrary value: maximum number of samples we require

min_nBAF_SNPs = 0

g_binwidth = 0.001
inputdir = "CN_joint_log2rs/"
bafdir = "BAF_joint_vals_fixed_purity/"
output_directory = "joint_analysis_fixed_purity/"
ploidyfile = "joint_ploidy.txt"
purityfile = "joint_purities.txt"
# read the filtered data that compares Xiaohong's segmentation data with raw SNP data

doubled = [[141, 21060], [141, 21062], [141, 21064], [163, 19208], [163, 19214], [194, 19868], [194, 19880], [450, 18974], [450, 18982], [512, 18744], [512, 18746], [512, 18748], [512, 18750], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944], [772, 18946], [848, 18794], [884, 20354], [884, 20358], [954, 20014], [954, 20016], [954, 20018], [991, 20600], [997, 20656], [997, 20666], [997, 20668], [997, 20672], [997, 20674], [1006, 21104], [1044, 20856], [1044, 20864], [997, 20658], [997, 20660], [660, 19266], [660, 19270], [740, 20000], [997, 20662], [997, 20664]]


def sort_log2r(outvec, call, intA, intB, all_data, double_loss_data, loss_data, wt_data, loh_data, single_gain_03_data, single_gain_12_data, balanced_gain_data, double_gain_04_data, double_gain_13_data, gain_data):
    all_data.append(outvec)
    if (call == 0):
        double_loss_data.append(outvec)
    elif (call == 1):
        loss_data.append(outvec)
    elif call == 2:
        if intA == intB:
            wt_data.append(outvec)
        else:
            loh_data.append(outvec)
    elif (call == 3):
        if abs(intA-intB)==3:
            single_gain_03_data.append(outvec)
        else:
            single_gain_12_data.append(outvec)
    elif call == 4:
        if intA==intB:
            balanced_gain_data.append(outvec)
        elif abs(intA-intB) == 4:
            double_gain_04_data.append(outvec)
        else:
            double_gain_13_data.append(outvec)
    else:
        gain_data.append(outvec)
    if (call > 2 and intA==intB):
        all_balanced_gain_data.append(outvec)

all_data = []
double_loss_data = []
loss_data = []
wt_data = []
loh_data = []
single_gain_03_data = []
single_gain_12_data = []
double_gain_04_data = []
double_gain_13_data = []
balanced_gain_data = []
all_balanced_gain_data = []
gain_data = []


ploidies = {}
pfile = open(ploidyfile, "r")
for line in pfile:
    if line.find("x") != -1:
        continue
    (id, ploidy) = line.rstrip().split()
    id = id[1:len(id)-1]
    (patient, sample) = id.split("_")
    ploidy = float(ploidy)
    #print "Setting", patient, sample, ploidy
    ploidies[(patient, sample)] = ploidy

purities = {}
pfile = open(purityfile, "r")
for line in pfile:
    if line.find("x") != -1:
        continue
    (id, purity) = line.rstrip().split()
    id = id[1:len(id)-1]
    (patient, sample) = id.split("_")
    purity = float(purity)
    #print "Setting", patient, sample, purity
    purities[(patient, sample)] = purity

for ps in ploidies:
    if ps in purities:
        print ps[0], ps[1], ploidies[ps], purities[ps]

flist = []
for (_, _, f) in walk(inputdir):
    flist += f
    
for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (int(patient), int(sample)) in doubled:
        print "Skipping due to doubling:", patient, sample
        continue
    ploidy = ploidies[(patient, sample)]
    purity = purities[(patient, sample)]
    print "Processing", patient, sample, ploidy, purity

    baffilename = patient + "_" + sample + "_avgbafvs.txt"
    baffile = open(bafdir + baffilename, "r")
    bafs = {}
    for line in baffile:
        (chr, start, end, rawA, rawB, intA, intB, avgBAF, nBAF_SNPs) = line.rstrip().split()
        if (chr=="chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            continue
        nBAF_SNPs = int(nBAF_SNPs)
        if (nBAF_SNPs < min_nBAF_SNPs):
            avgBAF = "NA"
        bafs[(chr, start, end)] = avgBAF

    total_n = 0
    log2rfile = open(inputdir + f, "r")
    for line in log2rfile:
        if line.find("chr") != -1:
            continue
        (chr, start, end, rawA, rawB, intA, intB, avg_log2r, nSNPs) = line.rstrip().split()
        if (chr == "chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            continue
        nSNPs = int(nSNPs)
        if use_min and nSNPs < nsamples_min:
            continue
        elif (use_max and nSNPs > nsamples_max):
            continue
        rawA = float(rawA)
        rawB = float(rawB)
        intA = int(intA)
        intB = int(intB)
        avg_log2r = float(avg_log2r)
        call = intA + intB
        bafdiff = bafs.get((chr, start, end))
        if (bafdiff == None):
            print "Can't find BAF difference for patient ", patient, " sample ", sample, ", segment ", chr, start, end
            continue
        if (bafdiff == "NA"):
            continue;
        bafdiff = float(bafdiff)
        disjoin = abs(rawA - intA) + abs(rawB - intB)
        outvec = [avg_log2r, "", "", "", "", "", ploidy, purity, nSNPs, disjoin]
        pindex = 5
        if purity>=0.95:
            pindex = 1
        elif purity >= 0.9:
            pindex = 2
        elif purity >= 0.7:
            pindex = 3
        elif purity >= 0.5:
            pindex = 4
        outvec[pindex] = bafdiff
        sort_log2r(outvec, call, intA, intB, all_data, double_loss_data, loss_data, wt_data, loh_data, single_gain_03_data, single_gain_12_data, balanced_gain_data, double_gain_04_data, double_gain_13_data, gain_data)
            

labels = "log2r\tBAF, 0.95<pur<=1\tBAF, 0.9<pur<0.95\tBAF, 0.7<pur<0.9\tBAF, 0.5<pur<0.7\tBAF, 0<pur<0.5\tploidy\tpurity\tnSNPs\tdisjoin"
lsl.saveScatterPlot(all_data, output_directory + "all_scatter.txt", labels)

lsl.saveScatterPlot(double_loss_data, output_directory + "double_loss_scatter.txt", labels)

lsl.saveScatterPlot(loss_data, output_directory + "loss_scatter.txt", labels)

lsl.saveScatterPlot(wt_data, output_directory + "wt_scatter.txt", labels)

lsl.saveScatterPlot(loh_data, output_directory + "loh_scatter.txt", labels)

lsl.saveScatterPlot(single_gain_03_data, output_directory + "single_gain_03_scatter.txt", labels)

lsl.saveScatterPlot(single_gain_12_data, output_directory + "single_gain_12_scatter.txt", labels)

lsl.saveScatterPlot(balanced_gain_data, output_directory + "balanced_gain_scatter.txt", labels)

lsl.saveScatterPlot(double_gain_04_data, output_directory + "double_gain_04_scatter.txt", labels)

lsl.saveScatterPlot(double_gain_13_data, output_directory + "double_gain_13_scatter.txt", labels)

lsl.saveScatterPlot(gain_data, output_directory + "other_gain_scatter.txt", labels)

lsl.saveScatterPlot(all_balanced_gain_data, output_directory +  "all_balanced_gain_scatter.txt", labels)

