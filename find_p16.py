# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:27:59 2016

@author: lpsmith
""" 

from __future__ import division
from os import walk
import lucianSNPLibrary as lsl

rejoin = ""
#comment or uncomment the following line to toggle whether the rejoined values are to be used.
rejoin = "_rejoined"

# read the processed data from the file Xiaohong sent me, which has segmentation and non-2N log2R values

#filenames = ["660_19272_avglog2rs.txt"]
filenames = []
for (_, _, f) in walk("CN_calc_log2rs" + rejoin + "/"):
    filenames += f
    break        

all_data = []
loss_data = []
double_loss_data = []
wt_data = []

outfilename = "short_segments/p16" + rejoin + ".txt"
outfile = open(outfilename, "w")
outline = "patient\tsample\tchr\tstart\tend\tXiao_log2r\tcall\tnum_log2rs\tavg(log2r)\tstdev\n"
outfile.write(outline)

for file in filenames:
    split= file.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    segmented_file = open("CN_calc_log2rs" + rejoin + "/" + file, "r")
        
    for line in segmented_file:
        if (line.find("chr") != -1):
            continue
        this_line = patient + "\t" + sample + "\t" + line
        (chr, start, end, xLog2r, call, nlog2r, log2r, stdev) = line.rstrip().split()
        chr = int(chr)
        if (chr != 9):
            continue
        start = int(start)
        if (start > 21995301):
            continue
        if (end=="inf"):
            end = 30000000 #Greater than the end of the gene
        end = int(end)
        if (end < 21967752):
            continue
        all_data.append(float(log2r))
        if call == "Loss":
            loss_data.append(float(log2r))
        elif call == "Double_d":
            double_loss_data.append(float(log2r))
        elif call == "wt":
            wt_data.append(float(log2r))
        outfile.write(this_line)

outfile.close()
lsl.createPrintAndSaveHistogram(double_loss_data, "short_segments/p16_double" + rejoin + ".txt", 0.001, axis=(-3.5, 1.5, 0))
lsl.createPrintAndSaveHistogram(loss_data, "short_segments/p16_loss" + rejoin + ".txt", 0.001, axis=(-3.5, 1.5, 0))
lsl.createPrintAndSaveHistogram(wt_data, "short_segments/p16_wt" + rejoin + ".txt", 0.001, axis=(-3.5, 1.5, 0))
lsl.createPrintAndSaveHistogram(all_data, "short_segments/p16_all" + rejoin + ".txt", 0.001, axis=(-3.5, 1.5, 0))
