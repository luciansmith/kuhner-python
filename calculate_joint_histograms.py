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

use_baf = False
baf_boundaries = [0.6, 0.7, 0.8, 0.9] #and an implied 'greater than that'
use_lengths = False
length_boundaries = [20, 100, 1000, 10000] #and an implied 'greater than that'

use_disjoin = True
disjoin_boundaries = [.2, .4, .6, .8]

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
bafdir = "BAF_joint_vals/"
output_directory = "joint_analysis/"
ploidyfile = "joint_ploidy.txt"
# read the filtered data that compares Xiaohong's segmentation data with raw SNP data

doubled = [[141, 21060], [141, 21062], [141, 21064], [163, 19208], [163, 19214], [194, 19868], [194, 19880], [450, 18974], [450, 18982], [512, 18744], [512, 18746], [512, 18748], [512, 18750], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944], [772, 18946], [848, 18794], [884, 20354], [884, 20358], [954, 20014], [954, 20016], [954, 20018], [991, 20600], [997, 20656], [997, 20666], [997, 20668], [997, 20672], [997, 20674], [1006, 21104], [1044, 20856], [1044, 20864], [997, 20658], [997, 20660], [660, 19266], [660, 19270], [740, 20000], [997, 20662], [997, 20664]]


def sort_log2r(avg_log2r, index, call, intA, intB, all_data, double_loss_data, loss_data, wt_data, loh_data, single_gain_03_data, single_gain_12_data, balanced_gain_data, double_gain_04_data, double_gain_13_data, gain_data):
    all_data[index].append(avg_log2r)
    all_data[5].append(avg_log2r)
    if (call == 0):
        double_loss_data[index].append(avg_log2r)
        double_loss_data[5].append(avg_log2r)
    elif (call == 1):
        loss_data[index].append(avg_log2r)
        loss_data[5].append(avg_log2r)
    elif call == 2:
        if intA == intB:
            wt_data[index].append(avg_log2r)
            wt_data[5].append(avg_log2r)
        else:
            loh_data[index].append(avg_log2r)
            loh_data[5].append(avg_log2r)
    elif (call == 3):
        if abs(intA-intB)==3:
            single_gain_03_data[index].append(avg_log2r)
            single_gain_03_data[5].append(avg_log2r)
        else:
            single_gain_12_data[index].append(avg_log2r)
            single_gain_12_data[5].append(avg_log2r)
    elif call == 4:
        if intA==intB:
            balanced_gain_data[index].append(avg_log2r)
            balanced_gain_data[5].append(avg_log2r)
        elif abs(intA-intB) == 4:
            double_gain_04_data[index].append(avg_log2r)
            double_gain_04_data[5].append(avg_log2r)
        else:
            double_gain_13_data[index].append(avg_log2r)
            double_gain_13_data[5].append(avg_log2r)
    else:
        gain_data[index].append(avg_log2r)
        gain_data[5].append(avg_log2r)
    if (call > 2 and intA==intB):
        all_balanced_gain_data[index].append(avg_log2r)
        all_balanced_gain_data[5].append(avg_log2r)

all_data = [[], []]
double_loss_data = [[], []]
loss_data = [[], []]
wt_data = [[], []]
loh_data = [[], []]
single_gain_data = [[], []]
double_gain_data = [[], []]
balanced_gain_data = [[], []]
all_balanced_gain_data = [[], []]
gain_data = [[], []]

if use_baf or use_lengths or use_disjoin:
    all_data = [[], [], [], [], [], []]
    double_loss_data = [[], [], [], [], [], []]
    loss_data = [[], [], [], [], [], []]
    wt_data = [[], [], [], [], [], []]
    loh_data = [[], [], [], [], [], []]
    single_gain_03_data = [[], [], [], [], [], []]
    single_gain_12_data = [[], [], [], [], [], []]
    double_gain_04_data = [[], [], [], [], [], []]
    double_gain_13_data = [[], [], [], [], [], []]
    balanced_gain_data = [[], [], [], [], [], []]
    all_balanced_gain_data = [[], [], [], [], [], []]
    gain_data = [[], [], [], [], [], []]


ploidies = {}
pfile = open(ploidyfile, "r")
for line in pfile:
    if line.find("x") != -1:
        continue
    (id, ploidy) = line.rstrip().split()
    id = id[1:len(id)-1]
    (patient, sample) = id.split("_")
    ploidy = float(ploidy)
    if (ploidy < 3):
        ploidy = 2
    else:
        ploidy = 4
    print "Setting", patient, sample, ploidy
    ploidies[(patient, sample)] = ploidy


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
    print "Processing", patient, sample, ploidy

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
        if (use_lengths):
            index = 4
            for n in range(len(length_boundaries)):
                if (nSNPs <= length_boundaries[n]):
                    index = n
                    break;
            sort_log2r(avg_log2r, index, call, intA, intB, all_data, double_loss_data, loss_data, wt_data, loh_data, single_gain_03_data, single_gain_12_data, balanced_gain_data, double_gain_04_data, double_gain_13_data, gain_data)
        elif (use_baf):
            index = 4
            for n in range(len(baf_boundaries)):
                if (bafdiff <= baf_boundaries[n]):
                    index = n
                    break;
            sort_log2r(avg_log2r, index, call, intA, intB, all_data, double_loss_data, loss_data, wt_data, loh_data, single_gain_03_data, single_gain_12_data, balanced_gain_data, double_gain_04_data, double_gain_13_data, gain_data)
        elif (use_disjoin):
            disjoin = abs(rawA - intA) + abs(rawB - intB)
            index = 4
            for n in range(len(disjoin_boundaries)):
                if (disjoin <= disjoin_boundaries[n]):
                    index = n
                    break;
            sort_log2r(avg_log2r, index, call, intA, intB, all_data, double_loss_data, loss_data, wt_data, loh_data, single_gain_03_data, single_gain_12_data, balanced_gain_data, double_gain_04_data, double_gain_13_data, gain_data)
        else:
            index=1
            if (ploidy == 2):
                index=0
            sort_log2r(avg_log2r, index, call, intA, intB, all_data, double_loss_data, loss_data, wt_data, loh_data, single_gain_03_data, single_gain_12_data, balanced_gain_data, double_gain_04_data, double_gain_13_data, gain_data)
            

rangestr = ""
ydata = []
if (use_baf):
    rangestr += "_bafs"
    ydata = ["0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0", "All BAFs"]
elif (use_lengths):
    rangestr += "_lengths"
    ydata = ["1-20", "21-100", "101-1000", "1001-10000", "10001+", "All lengths"]
elif (use_disjoin):
    rangestr += "_disjoins"
    ydata = disjoin_boundaries = ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8+", "All disjoins"]
else:
    rangestr += "_ploidies"
    ydata = ["ploidy_2", "ploidy_4"]

if (use_max):
    rangestr += "_only_" + str(nsamples_min) + "-" + str(nsamples_max)

boundaries = [-1, 1.2, 0]


print "Overall histograms:"
lsl.createPrintAndSaveMultipleHistograms(all_data, output_directory + "all_hist" + rangestr + ".txt", g_binwidth, ydata=ydata, axis=boundaries)

print "Double-loss histograms:"
lsl.createPrintAndSaveMultipleHistograms(double_loss_data, output_directory + "double_loss_hist" + rangestr + ".txt", g_binwidth, ydata=ydata)

print "Loss histograms:"
lsl.createPrintAndSaveMultipleHistograms(loss_data, output_directory + "loss_hist" + rangestr + ".txt", g_binwidth, ydata=ydata, axis=boundaries)

print "WT histograms:"
lsl.createPrintAndSaveMultipleHistograms(wt_data, output_directory + "wt_hist" + rangestr + ".txt", g_binwidth, ydata=ydata, axis=boundaries)

print "LOH histograms:"
lsl.createPrintAndSaveMultipleHistograms(loh_data, output_directory + "loh_hist" + rangestr + ".txt", g_binwidth, ydata=ydata, axis=boundaries)

print "Single Gain 0,3 histograms:"
lsl.createPrintAndSaveMultipleHistograms(single_gain_03_data, output_directory + "single_gain_03_hist" + rangestr + ".txt", g_binwidth, ydata=ydata,    axis=boundaries)

print "Single Gain 1,2 histograms:"
lsl.createPrintAndSaveMultipleHistograms(single_gain_12_data, output_directory + "single_gain_12_hist" + rangestr + ".txt", g_binwidth, ydata=ydata,    axis=boundaries)

print "Balanced gain histograms:"
lsl.createPrintAndSaveMultipleHistograms(balanced_gain_data, output_directory + "balanced_gain_hist" + rangestr + ".txt", g_binwidth, ydata=ydata,    axis=boundaries)

print "Double Gain 0,4 histograms:"
lsl.createPrintAndSaveMultipleHistograms(double_gain_04_data, output_directory + "double_gain_04_hist" + rangestr + ".txt", g_binwidth, ydata=ydata,    axis=boundaries)

print "Double Gain 1,3 histograms:"
lsl.createPrintAndSaveMultipleHistograms(double_gain_13_data, output_directory + "double_gain_13_hist" + rangestr + ".txt", g_binwidth, ydata=ydata,    axis=boundaries)

print "Other Gain histograms:"
lsl.createPrintAndSaveMultipleHistograms(gain_data, output_directory + "other_gain_hist" + rangestr + ".txt", g_binwidth, ydata=ydata, axis=boundaries)

print "All balanced gain histograms (A==B):"
lsl.createPrintAndSaveMultipleHistograms(all_balanced_gain_data, output_directory +  "all_balanced__gain_hist" + rangestr + ".txt", g_binwidth, ydata=ydata, axis=boundaries)

