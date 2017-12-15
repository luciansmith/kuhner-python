# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk

import lucianSNPLibrary as lsl
            
# read the filtered data that compares Xiaohong's segmentation data with raw SNP data

#filenames = ["1049_20780_avglog2rs.txt", "1049_20782_avglog2rs.txt"]
filenames = []
for (_, _, f) in walk("CN_calc_log2rs/"):
    filenames += f
    break        
doubled = [[141, 21060], [141, 21062], [141, 21064], [163, 19208], [163, 19214], [194, 19868], [194, 19880], [450, 18974], [450, 18982], [512, 18744], [512, 18746], [512, 18748], [512, 18750], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944], [772, 18946], [848, 18794], [884, 20354], [884, 20358], [954, 20014], [954, 20016], [954, 20018], [991, 20600], [997, 20656], [997, 20666], [997, 20668], [997, 20672], [997, 20674], [1006, 21104], [1044, 20856], [1044, 20864]]

new_doubled = [[868, 18714], [141, 21062], [146, 21358], [163, 19208], [163, 19214], [194, 19868], [509, 19000], [512, 18748], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944]]

#doubled += new_doubled

rejected_doubles = [[997, 20658], [997, 20660], [660, 19266], [660, 19270], [740, 20000], [997, 20662], [997, 20664]]

doubled += rejected_doubles
# All the doubled data was actually recalibrated.

combined_histograms = {}
combined_doubled_histograms = {}

double_loss_data = []
double_loss_from_doubled_data = []
loss_data = []
loss_from_doubled_data = []
wt_data = []
gain_data = []
balanced_gain_data = []
all_data = []

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
    file = open("CN_calc_log2rs/" + filename, "r")
    total_n = 0
    for line in file:
        (chr, start, end, x_log2r, call, n_log2r, avg_log2r, stdev) = line.rstrip().split()
        if (chr == "chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            continue
        if (start == "0"):
            continue
        if (end == "inf"):
            continue
        n_log2r = int(n_log2r)
        all_data.append(n_log2r)
        if (call == "Double_d"):
            if ([patient, sample] in doubled):
                double_loss_from_doubled_data.append(n_log2r)
            else:
                double_loss_data.append(n_log2r)
        elif (call == "Loss"):
            if ([patient, sample] in doubled):
                loss_from_doubled_data.append(n_log2r)
            else:
                loss_data.append(n_log2r)
            #loss_data.append(avg_log2r)
        elif (call == "wt"):
            wt_data.append(n_log2r)
        elif (call == "Gain"):
            gain_data.append(n_log2r)
        elif (call == "Balanced_gain"):
            balanced_gain_data.append(n_log2r)
        else:
            print "Unknown call ", call
    

binwidth = 0.001

print "Double-loss from doubled genomes histogram:"
lsl.createPrintAndSaveHistogram(double_loss_from_doubled_data, "CN_seglength_histograms/double_loss_from_doubled_hist.txt", binwidth)

print "Loss from doubled genomes histogram:"
lsl.createPrintAndSaveHistogram(loss_from_doubled_data, "CN_seglength_histograms/loss_from_doubled_hist.txt", binwidth)

print "Double-loss from doubled genomes histogram:"
lsl.createPrintAndSaveHistogram(double_loss_data, "CN_seglength_histograms/double_loss.txt", binwidth)

print "Loss from doubled genomes histogram:"
lsl.createPrintAndSaveHistogram(loss_data, "CN_seglength_histograms/loss.txt", binwidth)

print "WT histogram:"
lsl.createPrintAndSaveHistogram(wt_data, "CN_seglength_histograms/wt_hist.txt", "w")

print "Balanced gain histogram:"
lsl.createPrintAndSaveHistogram(balanced_gain_data,"CN_seglength_histograms/balanced_gain_hist.txt", "w")

print "Gain histogram:"
lsl.createPrintAndSaveHistogram(gain_data,"CN_seglength_histograms/gain_hist.txt", "w")

print "All data histogram:"
lsl.createPrintAndSaveHistogram(all_data,"CN_seglength_histograms/all_hist.txt", "w")
