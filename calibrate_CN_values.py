# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile
import numpy

nsamples_min = 5 #Arbitrary value: minimum number of samples we require

# read the filtered data that compares Xiaohong's segmentation data with raw SNP data

#filenames = ["1001_20802_avglog2rs.txt"]
filenames = []
for (_, _, f) in walk("CN_calc_log2rs/"):
    filenames += f
    break        

bgainvals = [[], [], [], [], []]
gainvals = [[], [], [], [], []]
wtvals = [[], [], [], [], []]
lossvals = [[], [], [], [], []]
ddelvals = [[], [], [], [], []]
avg_err = []
lavg_err = []
avgsh_err = []
lavgsh_err = []
for filename in filenames:
    split= filename.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    file = open("CN_calc_log2rs/" + filename, "r")
    for line in file:
        (chr, start, end, x_log2r, call, n_log2r, fiveperc, avg_log2r, lavg_log2r, avg_log2r_short, lavg_log2r_short) = line.rstrip().split()
        if (chr == "chr"):
            continue
        n_log2r = int(n_log2r)
        if (n_log2r < nsamples_min):
            continue
        x_log2r = float(x_log2r)
        avg_log2r = float(avg_log2r)
        avg_log2r_short = float(avg_log2r_short)
        lavg_log2r = float(lavg_log2r)
        lavg_log2r_short = float(lavg_log2r_short)
        if (call != "wt"):
            avg_err.append(abs(x_log2r - avg_log2r))
            lavg_err.append(abs(x_log2r - lavg_log2r))
            avgsh_err.append(abs(x_log2r - avg_log2r_short))
            lavgsh_err.append(abs(x_log2r - lavg_log2r_short))
        if (call == "Balanced_gain"):
            bgainvals[0].append(x_log2r)
            bgainvals[1].append(avg_log2r)
            bgainvals[2].append(lavg_log2r)
            bgainvals[3].append(avg_log2r_short)
            bgainvals[4].append(lavg_log2r_short)
        elif (call == "Gain"):
            gainvals[0].append(x_log2r)
            gainvals[1].append(avg_log2r)
            gainvals[2].append(lavg_log2r)
            gainvals[3].append(avg_log2r_short)
            gainvals[4].append(lavg_log2r_short)
        elif (call == "wt"):
            wtvals[0].append(x_log2r)
            wtvals[1].append(avg_log2r)
            wtvals[2].append(lavg_log2r)
            wtvals[3].append(avg_log2r_short)
            wtvals[4].append(lavg_log2r_short)
        elif (call == "Loss"):
            lossvals[0].append(x_log2r)
            lossvals[1].append(avg_log2r)
            lossvals[2].append(lavg_log2r)
            lossvals[3].append(avg_log2r_short)
            lossvals[4].append(lavg_log2r_short)
        elif (call == "Double_d"):
            ddelvals[0].append(x_log2r)
            ddelvals[1].append(avg_log2r)
            ddelvals[2].append(lavg_log2r)
            ddelvals[3].append(avg_log2r_short)
            ddelvals[4].append(lavg_log2r_short)
        else:
            print "unknown call '" + call + "'."
        
print "Error between Xiaohong call and straight log2r average: ", str(numpy.average(avg_err))
print "Error between Xiaohong call and truncated log2r average: ", str(numpy.average(avgsh_err))
print "Error between Xiaohong call and average of 2^log2r's: ", str(numpy.average(lavg_err))
print "Error between Xiaohong call and truncated average of 2^log2r's: ", str(numpy.average(lavgsh_err))

print ""
print "Balanced gain log2r averages: ", str(numpy.average(bgainvals[1])), ", ", str(numpy.average(bgainvals[2])), ", ", str(numpy.average(bgainvals[3])), ", ", str(numpy.average(bgainvals[4]))
print "Gain log2r averages: ", str(numpy.average(gainvals[1])), ", ", str(numpy.average(gainvals[2])), ", ", str(numpy.average(gainvals[3])), ", ", str(numpy.average(gainvals[4]))
print "wt log2r averages: ", str(numpy.average(wtvals[1])), ", ", str(numpy.average(wtvals[2])), ", ", str(numpy.average(wtvals[3])), ", ", str(numpy.average(wtvals[4]))
print "Loss log2r averages: ", str(numpy.average(lossvals[1])), ", ", str(numpy.average(lossvals[2])), ", ", str(numpy.average(lossvals[3])), ", ", str(numpy.average(lossvals[4]))
print "Double deletion log2r averages: ", str(numpy.average(ddelvals[1])), ", ", str(numpy.average(ddelvals[2])), ", ", str(numpy.average(ddelvals[3])), ", ", str(numpy.average(ddelvals[4]))
print ""
print "Balanced gain log2r medians: ", str(numpy.median(bgainvals[1])), ", ", str(numpy.median(bgainvals[2])), ", ", str(numpy.median(bgainvals[3])), ", ", str(numpy.median(bgainvals[4]))
print "Gain log2r medians: ", str(numpy.median(gainvals[1])), ", ", str(numpy.median(gainvals[2])), ", ", str(numpy.median(gainvals[3])), ", ", str(numpy.median(gainvals[4]))
print "wt log2r medians: ", str(numpy.median(wtvals[1])), ", ", str(numpy.median(wtvals[2])), ", ", str(numpy.median(wtvals[3])), ", ", str(numpy.median(wtvals[4]))
print "Loss log2r medians: ", str(numpy.median(lossvals[1])), ", ", str(numpy.median(lossvals[2])), ", ", str(numpy.median(lossvals[3])), ", ", str(numpy.median(lossvals[4]))
print "Double deletion log2r medians: ", str(numpy.median(ddelvals[1])), ", ", str(numpy.median(ddelvals[2])), ", ", str(numpy.median(ddelvals[3])), ", ", str(numpy.median(ddelvals[4]))
print ""
print "Balanced gain log2r max: ", str(max(bgainvals[1])), ", ", str(max(bgainvals[2])), ", ", str(max(bgainvals[3])), ", ", str(max(bgainvals[4]))
print "Gain log2r max: ", str(max(gainvals[1])), ", ", str(max(gainvals[2])), ", ", str(max(gainvals[3])), ", ", str(max(gainvals[4]))
print "wt log2r max: ", str(max(wtvals[1])), ", ", str(max(wtvals[2])), ", ", str(max(wtvals[3])), ", ", str(max(wtvals[4]))
print "Loss log2r max: ", str(max(lossvals[1])), ", ", str(max(lossvals[2])), ", ", str(max(lossvals[3])), ", ", str(max(lossvals[4]))
print "Double deletion log2r max: ", str(max(ddelvals[1])), ", ", str(max(ddelvals[2])), ", ", str(max(ddelvals[3])), ", ", str(max(ddelvals[4]))
print ""
print "Balanced gain log2r min: ", str(min(bgainvals[1])), ", ", str(min(bgainvals[2])), ", ", str(min(bgainvals[3])), ", ", str(min(bgainvals[4]))
print "Gain log2r min: ", str(min(gainvals[1])), ", ", str(min(gainvals[2])), ", ", str(min(gainvals[3])), ", ", str(min(gainvals[4]))
print "wt log2r min: ", str(min(wtvals[1])), ", ", str(min(wtvals[2])), ", ", str(min(wtvals[3])), ", ", str(min(wtvals[4]))
print "Loss log2r min: ", str(min(lossvals[1])), ", ", str(min(lossvals[2])), ", ", str(min(lossvals[3])), ", ", str(min(lossvals[4]))
print "Double deletion log2r min: ", str(min(ddelvals[1])), ", ", str(min(ddelvals[2])), ", ", str(min(ddelvals[3])), ", ", str(min(ddelvals[4]))

bgOut = open("Calibration_balanced_gain.txt", "w")
gOut  = open("Calibration_gain.txt", "w")
wtOut = open("Calibration_wt.txt", "w")
lOut  = open("Calibration_loss.txt", "w")
ddOut = open("Calibration_del.txt", "w")
outline = "average\taverage_short\t,logavg\tlogavg_short\n"
bgOut.write(outline)
gOut.write(outline)
wtOut.write(outline)
lOut.write(outline)
ddOut.write(outline)

for i in range(0, len(bgainvals[0])):
    bgOut.write(str(bgainvals[1][i]) + "\t" + str(bgainvals[2][i]) + "\t" + str(bgainvals[3][i])+ "\t" + str(bgainvals[4][i]) + "\n")
for i in range(0, len(gainvals[0])):
    gOut.write(str(gainvals[1][i]) + "\t" + str(gainvals[2][i]) + "\t" + str(gainvals[3][i]) + "\t" + str(gainvals[4][i]) + "\n")
for i in range(0, len(wtvals[0])):
    wtOut.write(str(wtvals[1][i]) + "\t" + str(wtvals[2][i]) + "\t" + str(wtvals[3][i]) +  "\t" + str(wtvals[4][i]) + "\n")
for i in range(0, len(lossvals[0])):
    lOut.write(str(lossvals[1][i]) + "\t" + str(lossvals[2][i]) + "\t" + str(lossvals[3][i]) +  "\t" + str(lossvals[4][i]) + "\n")
for i in range(0, len(ddelvals[0])):
    ddOut.write(str(ddelvals[1][i]) + "\t" + str(ddelvals[2][i]) + "\t" + str(ddelvals[3][i]) +  "\t" + str(ddelvals[4][i]) + "\n")
    