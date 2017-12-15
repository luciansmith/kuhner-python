# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division

import lucianSNPLibrary as lsl

from os import walk

nsamples_min = 21 #Arbitrary value: minimum number of samples we require
nsamples_max = 100000000

binwidth = 0.001
#indir = "CN_calc_log2rs/"
#outdir = "CN_smoothed_histograms/"
indir = "CN_calc_log2rs_rejoined/"
outdir = "CN_rejoined_histograms/"
srange = ""
if (nsamples_max > 0):
    srange = "_" + str(nsamples_min) + "-" + str(nsamples_max)

# read the filtered data that compares Xiaohong's segmentation data with raw SNP data

#filenames = ["954_20016_avglog2rs.txt", "1049_20782_avglog2rs.txt"]
filenames = []
for (_, _, f) in walk(indir):
    filenames += f
    break        

doubled = [[141, 21060], [141, 21062], [141, 21064], [163, 19208], [163, 19214], [194, 19868], [194, 19880], [450, 18974], [450, 18982], [512, 18744], [512, 18746], [512, 18748], [512, 18750], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944], [772, 18946], [848, 18794], [884, 20354], [884, 20358], [954, 20014], [954, 20016], [954, 20018], [991, 20600], [997, 20656], [997, 20666], [997, 20668], [997, 20672], [997, 20674], [1006, 21104], [1044, 20856], [1044, 20864], [997, 20658], [997, 20660], [997, 20662], [660, 19266], [660, 19270], [997, 20664], [740, 20000]]

#new_doubled = [[141, 21062], [163, 19208], [163, 19214], [194, 19868], [509, 19000], [512, 18748], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944]]

#doubled += new_doubled

#rejected_doubles = []

#doubled += rejected_doubles

all_data = []
double_loss_data = []
double_loss_from_doubled_data = []
loss_data = []
loss_from_doubled_data = []
wt_data = []
gain_data = []
balanced_gain_data = []

for filename in filenames:
    if (filename.find(".txt") == -1):
        continue
    split= filename.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    patient = int(patient)
    sample = int(sample)
    file = open(indir + filename, "r")
    total_n = 0
    sample_data = []
    for line in file:
        (chr, start, end, x_log2r, call, n_log2r, avg_log2r, stdev) = line.rstrip().split()
        if (chr == "chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            continue
        n_log2r = int(n_log2r)
        if (n_log2r < nsamples_min):
            continue
        if (nsamples_max > 0 and n_log2r > nsamples_max):
            continue
        total_n += n_log2r
        avg_log2r = float(avg_log2r)
        all_data.append(avg_log2r)
        sample_data.append(avg_log2r)
        if (call == "Double_d"):
            if ([patient, sample] in doubled):
                double_loss_from_doubled_data.append(avg_log2r)
            else:
                double_loss_data.append(avg_log2r)
        elif (call == "Loss"):
            if ([patient, sample] in doubled):
                loss_from_doubled_data.append(avg_log2r)
            else:
                loss_data.append(avg_log2r)
            #loss_data.append(avg_log2r)
        elif (call == "wt"):
            wt_data.append(avg_log2r)
        elif (call == "Gain"):
            gain_data.append(avg_log2r)
        elif (call == "Balanced_gain"):
            balanced_gain_data.append(avg_log2r)
        else:
            print "Unknown call ", call

    lsl.createPrintAndSaveHistogram(double_loss_from_doubled_data, outdir + str(patient) + "_" + str(sample) + "_smoothhist.txt", binwidth, show=False)


print "Double-loss from doubled genomes histogram:"
lsl.createPrintAndSaveHistogram(double_loss_from_doubled_data, outdir + "double_loss_from_doubled_hist" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))

print "Loss from doubled genomes histogram:"
lsl.createPrintAndSaveHistogram(loss_from_doubled_data, outdir + "loss_from_doubled_hist" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))

print "Double-loss histogram:"
lsl.createPrintAndSaveHistogram(double_loss_data, outdir + "double_loss" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))

print "Loss histogram:"
lsl.createPrintAndSaveHistogram(loss_data, outdir + "loss" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))

print "WT histogram:"
lsl.createPrintAndSaveHistogram(wt_data, outdir + "wt_hist" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))

print "Balanced gain histogram:"
lsl.createPrintAndSaveHistogram(balanced_gain_data,outdir + "balanced_gain_hist" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))

print "Gain histogram:"
lsl.createPrintAndSaveHistogram(gain_data,outdir + "gain_hist" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))

print "All data histogram:"
lsl.createPrintAndSaveHistogram(all_data,outdir + "all_hist" + srange + ".txt", binwidth, axis=(-3.5, 1.5, 0))
