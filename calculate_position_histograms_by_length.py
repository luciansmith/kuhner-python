# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk

import lucianSNPLibrary as lsl


nsamples_min = 10 #Arbitrary value: minimum number of samples we require
#nsamples_max = 12

# read the filtered data that compares Xiaohong's segmentation data with raw SNP data

#filenames = ["954_20016_avglog2rs.txt", "1049_20782_avglog2rs.txt"]
filenames = []
for (_, _, f) in walk("CN_calc_log2rs/"):
    filenames += f
    break        

all_data = []
len10_data = []
len13_data = []
len21_data = []
len51_data = []
len501_data = []
len5001_data = []
len50001_data = []

all_chr_data = []
len10_chr_data = []
len13_chr_data = []
len21_chr_data = []
len51_chr_data = []
len501_chr_data = []
len5001_chr_data = []
len50001_chr_data = []

for c in range(0,22):
    all_chr_data.append([])
    len10_chr_data.append([])
    len13_chr_data.append([])
    len21_chr_data.append([])
    len51_chr_data.append([])
    len501_chr_data.append([])
    len5001_chr_data.append([])
    len50001_chr_data.append([])
    

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
    for line in file:
        (chr, start, end, x_log2r, call, n_log2r, fiveperc, avg_log2r, lavg_log2r, avg_log2r_short, lavg_log2r_short) = line.rstrip().split()
        if (chr == "chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            continue
        n_log2r = int(n_log2r)
        if (n_log2r < nsamples_min):
            continue
        
        pos = lsl.getLengthFrom(chr, start, end)
        pos = int(start) + pos/2
        gpos = lsl.getGlobalPosition(chr, start, end)
        if (gpos<0):
            print "Got position", gpos, "from start points",chr, start, end
        all_data.append(gpos)
        n_log2r = int(n_log2r)
        if (n_log2r >= 10 and n_log2r <= 12):
            len10_data.append(gpos)
            len10_chr_data[chr-1].append(pos)
        elif (n_log2r >= 13 and n_log2r <= 20):
            len13_data.append(gpos)
            len13_chr_data[chr-1].append(pos)
        elif (n_log2r >= 21 and n_log2r <= 50):
            len21_data.append(gpos)
            len21_chr_data[chr-1].append(pos)
        elif (n_log2r >= 51 and n_log2r <= 500):
            len51_data.append(gpos)
            len51_chr_data[chr-1].append(pos)
        elif (n_log2r >= 501 and n_log2r <= 5000):
            len501_data.append(gpos)
            len501_chr_data[chr-1].append(pos)
        elif (n_log2r >= 5001 and n_log2r <= 50000):
            len5001_data.append(gpos)
            len5001_chr_data[chr-1].append(pos)
        elif (n_log2r >= 50001):
            len50001_data.append(gpos)
            len50001_chr_data[chr-1].append(pos)
        else:
            print "Error: uncaught n_log2r value", n_log2r
        
        
directory = "seglength_chr_histograms/"
binwidth = 1000000
lsl.createPrintAndSaveHistogram(len10_data, directory + "len10_data.txt", binwidth, xdata="Position")
lsl.createPrintAndSaveHistogram(len13_data, directory + "len13_data.txt", binwidth, xdata="Position")
lsl.createPrintAndSaveHistogram(len21_data, directory + "len21_data.txt", binwidth, xdata="Position")
lsl.createPrintAndSaveHistogram(len51_data, directory + "len51_data.txt", binwidth, xdata="Position")
lsl.createPrintAndSaveHistogram(len501_data, directory + "len501_data.txt", binwidth, xdata="Position")
lsl.createPrintAndSaveHistogram(len5001_data, directory + "len5001_data.txt", binwidth, xdata="Position")
lsl.createPrintAndSaveHistogram(len50001_data, directory + "len50001_data.txt", binwidth, xdata="Position")
lsl.createPrintAndSaveHistogram(all_data, directory + "all_data.txt", binwidth, xdata="Position")

for c in range(0,22):
    lsl.createPrintAndSaveHistogram(len10_chr_data[c], directory + "len10_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
    lsl.createPrintAndSaveHistogram(len13_chr_data[c], directory + "len13_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
    lsl.createPrintAndSaveHistogram(len21_chr_data[c], directory + "len21_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
    lsl.createPrintAndSaveHistogram(len51_chr_data[c], directory + "len51_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
    lsl.createPrintAndSaveHistogram(len501_chr_data[c], directory + "len501_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
    lsl.createPrintAndSaveHistogram(len5001_chr_data[c], directory + "len5001_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
    lsl.createPrintAndSaveHistogram(len50001_chr_data[c], directory + "len50001_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
    lsl.createPrintAndSaveHistogram(all_chr_data[c], directory + "all_chr_data_" + str(c+1) + ".txt", binwidth, xdata="Position")
