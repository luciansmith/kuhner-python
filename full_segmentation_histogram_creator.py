# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:43:28 2016

@author: lpsmith
"""

from __future__ import division
from os import walk

import lucianSNPLibrary as lsl

nsamples_min = 10 #Arbitrary value: minimum number of samples we require
data10_12 = []
data13_20 = []
data21_50 = []
data51_500 = []
data501_5000 = []
data5001_50000 = []
data50001_plus = []
dataall =[]

#fullseg_filenames = ["three_formal_cy_omni_mix3_b37RB.txt"]
fullseg_filenames = []
for (_, _, f) in walk("full_segmentation_output/"):
    fullseg_filenames += f
    break

discrepancies = open("full_segmentation_histograms/discrepancies.txt", "w")
for file in fullseg_filenames:
    handle = open("full_segmentation_output/" + file, "r")
    for line in handle:
        (chr, start, end, pmean, pnmarkers, nmarkers, meanlog2r) = line.rstrip().split("\t")
        if (chr=="chr"):
            continue
        if (pnmarkers != "?"):
            pnmarkers = int(pnmarkers)
            nmarkers = int(nmarkers)
            if (pnmarkers != nmarkers):
                print "Anomaly in", file, ": different nmarkers from partek vs. raw SNP data:"
                print "  ", line
                line = file + "\t" + line
                discrepancies.write(line)
        if (nmarkers < nsamples_min):
            continue
        meanlog2r = float(meanlog2r)
        dataall.append(meanlog2r)
        if (nmarkers < 13):
            data10_12.append(meanlog2r)
        elif (nmarkers < 21):
            data13_20.append(meanlog2r)
        elif (nmarkers < 51):
            data21_50.append(meanlog2r)
        elif (nmarkers < 501):
            data51_500.append(meanlog2r)
        elif (nmarkers < 5001):
            data501_5000.append(meanlog2r)
        elif (nmarkers < 50001):
            data5001_50000.append(meanlog2r)
        elif (nmarkers < 500001):
            data50001_plus.append(meanlog2r)

binwidth = 0.001
lsl.createPrintAndSaveHistogram(data10_12, "full_segmentation_histograms/data10_12.txt", binwidth)
lsl.createPrintAndSaveHistogram(data13_20, "full_segmentation_histograms/data13_20.txt", binwidth)
lsl.createPrintAndSaveHistogram(data21_50, "full_segmentation_histograms/data21_50.txt", binwidth)
lsl.createPrintAndSaveHistogram(data51_500, "full_segmentation_histograms/data51_500.txt", binwidth)
lsl.createPrintAndSaveHistogram(data501_5000, "full_segmentation_histograms/data501_5000.txt", binwidth)
lsl.createPrintAndSaveHistogram(data5001_50000, "full_segmentation_histograms/data5001_50000.txt", binwidth)
lsl.createPrintAndSaveHistogram(data50001_plus, "full_segmentation_histograms/data50001_plus.txt", binwidth)
lsl.createPrintAndSaveHistogram(dataall, "full_segmentation_histograms/dataall.txt", binwidth)
