# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk

import lucianSNPLibrary as lsl


use_baf = True

use_length = True
length_min = 1000000001
length_max = 10000000000

nsamples_min = 10 #Arbitrary value: minimum number of samples we require
use_max = False
nsamples_max = 500 #Arbitrary value: minimum number of samples we require
use_baf = False

g_binwidth = 0.001
output_directory = "presentation2/"

# read the filtered data that compares Xiaohong's segmentation data with raw SNP data


flist = []
SNPfiles = []
#SNPfiles.append(["954", "20014", "954_20014_segment_BAFdiff.txt"])
for (_, _, f) in walk("SNP_segment_BAFs/"):
    flist += f
   
for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (sample == "blood"):
        continue
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    SNPfiles += [[patient, sample, f]]

overlap_files = []
flist = []
for (_, _, f) in walk("CN_calc_log2rs/"):
    flist += f
    
for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    for snpfile in SNPfiles:
        if (patient == snpfile[0] and sample == snpfile[1]):
            overlap_files += [[f, snpfile[2]]]
            break
        

doubled = [[141, 21060], [141, 21062], [141, 21064], [163, 19208], [163, 19214], [194, 19868], [194, 19880], [450, 18974], [450, 18982], [512, 18744], [512, 18746], [512, 18748], [512, 18750], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944], [772, 18946], [848, 18794], [884, 20354], [884, 20358], [954, 20014], [954, 20016], [954, 20018], [991, 20600], [997, 20656], [997, 20666], [997, 20668], [997, 20672], [997, 20674], [1006, 21104], [1044, 20856], [1044, 20864]]

new_doubled = [[868, 18714], [141, 21062], [146, 21358], [163, 19208], [163, 19214], [194, 19868], [509, 19000], [512, 18748], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944]]

#doubled += new_doubled

rejected_doubles = [[997, 20658], [997, 20660], [660, 19266], [660, 19270], [740, 20000], [997, 20662], [997, 20664]]

doubled += rejected_doubles #Put these in anyway--they come from recalibrated data

double_loss_data = []
loss_data = []
wt_data = []
gain_data = []
balanced_gain_data = []
all_data = []
position_data= []

if use_baf:
    double_loss_data = [[], [], [], [], []]
    loss_data = [[], [], [], [], []]
    wt_data = [[], [], [], [], []]
    gain_data = [[], [], [], [], []]
    balanced_gain_data = [[], [], [], [], []]
    all_data = [[], [], [], [], []]

maxes = []

for filepair in overlap_files:
    log2rfilename = filepair[0]
    baffilename = filepair[1]
    alldata = {}
    split = log2rfilename.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    patient = int(patient)
    sample = int(sample)
    baffile = open("SNP_segment_BAFs/" + baffilename, "r")
    bafs = {}
    for line in baffile:
        (chr, start, end, avg, maximum, minimum, number) = line.rstrip().split()
        if (chr=="chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            continue
        bafs[(chr, start, end)] = avg

    total_n = 0
    log2rfile = open("CN_calc_log2rs/" + log2rfilename, "r")
    for line in log2rfile:
        (chr, start, end, x_log2r, call, n_log2r, fiveperc, avg_log2r, lavg_log2r, avg_log2r_short, lavg_log2r_short) = line.rstrip().split()
        if (chr == "chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            continue
        n_log2r = int(n_log2r)
        if (n_log2r < nsamples_min):
            continue
        if (use_length):
            length = lsl.getLengthFrom(chr, start, end)
            if (length < length_min):
                continue
            if (length > length_max):
                continue
        elif (use_max and n_log2r > nsamples_max):
                continue
        avg_log2r = float(avg_log2r)
        if (use_baf):
            bafdiff = bafs.get((chr, start, end))
            if (bafdiff == None):
                print "Can't find BAF difference for patient ", patient, " sample ", sample, ", segment ", chr, start, end
                continue
            if (bafdiff == "--"):
                continue;
            bafdiff = float(bafdiff)
            index = 4
            if (bafdiff <= 0.1):
                index = 0
            elif (bafdiff <= 0.2):
                index = 1
            elif (bafdiff <= 0.3):
                index = 2
            elif (bafdiff <= 0.4):
                index = 3

            all_data[index].append(avg_log2r)
            if (call == "Double_d"):
                double_loss_data[index].append(avg_log2r)
            elif (call == "Loss"):
                loss_data[index].append(avg_log2r)
            elif (call == "wt"):
                wt_data[index].append(avg_log2r)
            elif (call == "Gain"):
                gain_data[index].append(avg_log2r)
            elif (call == "Balanced_gain"):
                balanced_gain_data[index].append(avg_log2r)
            else:
                print "Unknown call ", call
        else:
            all_data.append(avg_log2r)
            if (call == "Double_d"):
                double_loss_data.append(avg_log2r)
            elif (call == "Loss"):
                loss_data.append(avg_log2r)
            elif (call == "wt"):
                wt_data.append(avg_log2r)
            elif (call == "Gain"):
                gain_data.append(avg_log2r)
            elif (call == "Balanced_gain"):
                balanced_gain_data.append(avg_log2r)
            else:
                print "Unknown call ", call


rangestr = "_"
if (use_max):
    rangestr += "only_" + str(nsamples_min) + "-" + str(nsamples_max) + "_"
if (use_length):
    rangestr = "_only_" + str(length_min) + "-" + str(length_max) + "_"

if (use_baf):
    print "Double-loss histograms:"
    index = 0
    combined_data = []
    for dataset in double_loss_data:
        lsl.createPrintAndSaveHistogram(dataset, output_directory + "double_loss_hist" + rangestr + str(index) + ".txt", g_binwidth)
        combined_data += dataset
        index += 1
    lsl.createPrintAndSaveHistogram(combined_data, output_directory + "double_loss_hist" + rangestr + "all.txt", g_binwidth)
    
    combined_data = []
    print "Loss histograms:"
    index = 0
    for dataset in loss_data:
        lsl.createPrintAndSaveHistogram(dataset, output_directory + "loss_hist" + rangestr + str(index) + ".txt", g_binwidth)
        combined_data += dataset
        index += 1
    lsl.createPrintAndSaveHistogram(combined_data, output_directory + "loss_hist" + rangestr + "all.txt", g_binwidth)
    
    combined_data = []
    print "WT histograms:"
    index = 0
    for dataset in wt_data:
        lsl.createPrintAndSaveHistogram(dataset, output_directory + "wt_hist" + rangestr + str(index) + ".txt", g_binwidth)
        combined_data += dataset
        index += 1
    lsl.createPrintAndSaveHistogram(combined_data, output_directory + "wt_hist" + rangestr + "all.txt", g_binwidth)
    
    combined_data = []
    print "Gain histograms:"
    index = 0
    for dataset in gain_data:
        lsl.createPrintAndSaveHistogram(dataset, output_directory + "gain_hist" + rangestr + str(index) + ".txt", g_binwidth)
        combined_data += dataset
        index += 1
    lsl.createPrintAndSaveHistogram(combined_data, output_directory + "gain_hist" + rangestr + "all.txt", g_binwidth)
    
    combined_data = []
    print "Balanced gain histograms:"
    index = 0
    for dataset in balanced_gain_data:
        lsl.createPrintAndSaveHistogram(dataset, output_directory + "balanced_gain_hist" + rangestr + str(index) + ".txt", g_binwidth)
        combined_data += dataset
        index += 1
    lsl.createPrintAndSaveHistogram(combined_data, output_directory + "balanced_gain_hist" + rangestr + "all.txt", g_binwidth)
    
    combined_data = []
    print "Overall histograms:"
    index = 0
    for dataset in all_data:
        lsl.createPrintAndSaveHistogram(dataset, output_directory + "all_hist" + rangestr + str(index) + ".txt", g_binwidth)
        combined_data += dataset
        index += 1
    #combined_data += leftover_data
    lsl.createPrintAndSaveHistogram(combined_data, output_directory + "all_hist" + rangestr + "all.txt", g_binwidth)
    print "Total number of points for range", rangestr, ": ", len(combined_data)

else:
    print "Double-loss histograms:"
    lsl.createPrintAndSaveHistogram(double_loss_data, output_directory + "double_loss_hist" + rangestr + ".txt", g_binwidth)
    
    print "Loss histograms:"
    lsl.createPrintAndSaveHistogram(loss_data, output_directory + "loss_hist" + rangestr + ".txt", g_binwidth)
    
    print "WT histograms:"
    lsl.createPrintAndSaveHistogram(wt_data, output_directory + "wt_hist" + rangestr + ".txt", g_binwidth)
    
    print "Gain histograms:"
    lsl.createPrintAndSaveHistogram(gain_data, output_directory + "gain_hist" + rangestr + ".txt", g_binwidth)
    
    print "Balanced gain histograms:"
    lsl.createPrintAndSaveHistogram(balanced_gain_data, output_directory + "balanced_gain_hist" + rangestr + ".txt", g_binwidth)
    
    print "Overall histograms:"
    lsl.createPrintAndSaveHistogram(all_data, output_directory + "all_hist" + rangestr + ".txt", g_binwidth)
    

