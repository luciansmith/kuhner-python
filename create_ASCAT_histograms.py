# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
#from os.path import isfile
import lucianSNPLibrary as lsl

g_binwidth = 0.001
nsamples_min = 0 #Arbitrary value: minimum number of samples we require
nsamples_max = 9000000000

doubled = {(141, 21060), (141, 21062), (141, 21064), (163, 19208), (163, 19214), (194, 19868), (194, 19880), (450, 18974), (450, 18982), (512, 18744), (512, 18746), (512, 18748), (512, 18750), (512, 18762), (660, 19260), (660, 19262), (660, 19264), (660, 19272), (664, 19954), (772, 18944), (772, 18946), (848, 18794), (884, 20354), (884, 20358), (954, 20014), (954, 20016), (954, 20018), (991, 20600), (997, 20656), (997, 20666), (997, 20668), (997, 20672), (997, 20674), (1006, 21104), (1044, 20856), (1044, 20864)}

#new_doubled = [[868, 18714], [141, 21062], [146, 21358], [163, 19208], [163, 19214], [194, 19868], [509, 19000], [512, 18748], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944]]

#doubled += new_doubled

#rejected_doubles = [[997, 20658], [997, 20660], [660, 19266], [660, 19270], [740, 20000], [997, 20662], [997, 20664]]

#doubled += rejected_doubles

all_data = []
double_loss_data = []
loss_data = []
wt_data = []
LOH_data = []
gain_data = []
gain3_data = []
gain4_data = []
gain5_data = []
gain6_data = []

data_00 = []
data_01 = []
data_02 = []
data_03 = []
data_11 = []
data_12 = []
data_13 = []
data_22 = []
data_23 = []
data_33 = []

ascatfile = open("all_lesions.txt", "r")
for line in ascatfile:
    (patient, sample, chr, start, end, nsnps, log2r, intA, intB, floatA, floatB) = line.rstrip().split()
    if (chr == "chrom"):
        continue
    chr = int(chr)
    if (chr >= 23):
        continue
    nsnps = int(nsnps)
    if (nsnps < nsamples_min):
        continue
    if (nsnps > nsamples_max):
        continue
    if (log2r == "None"):
        continue
    if (int(patient), int(sample)) in doubled:
        #print "Skipping", patient, sample
        continue
    log2r = float(log2r)
    intA = int(intA)
    intB = int(intB)
    floatA = float(floatA)
    floatB = float(floatB)
        
    all_data.append(log2r)
    CN = intA + intB
    if (CN == 0):
        double_loss_data.append(log2r)
    elif (CN == 1):
        loss_data.append(log2r)
    elif (CN == 2):
        wt_data.append(log2r)
        if (intA == 2 or intB == 2):
            LOH_data.append(log2r)
    else:
        gain_data.append(log2r)
        if CN == 3:
            gain3_data.append(log2r)
        elif CN == 4:
            gain4_data.append(log2r)
        elif CN == 5:
            gain5_data.append(log2r)
        elif CN == 6:
            gain6_data.append(log2r)
    if intB == 0:
        if intA==0:
            data_00.append(log2r)
        elif intA==1:
            data_01.append(log2r)
        elif intA==2:
            data_02.append(log2r)
        elif intA==3:
            data_03.append(log2r)
    elif intB == 1:
        if intA==1:
            data_11.append(log2r)
        elif intA==2:
            data_12.append(log2r)
        elif intA==3:
            data_13.append(log2r)
    elif intB == 2:
        if intA==2:
            data_22.append(log2r)
        elif intA==3:
            data_23.append(log2r)
    elif intB == 3:
        if intA==3:
            data_33.append(log2r)

rangestr = "_" + str(nsamples_min) + "-" + str(nsamples_max) + "_"

thisaxis=[-3.5, 1.5, 0]
                    
print "Double-loss histograms:"
lsl.createPrintAndSaveHistogram(double_loss_data, "ASCAT_smoothed_histograms/double_loss_hist_a" + rangestr + ".txt", g_binwidth, axis=[-3.5, 1.5, 0])

print "Loss histograms:"
lsl.createPrintAndSaveHistogram(loss_data, "ASCAT_smoothed_histograms/loss_hist_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "WT histograms:"
lsl.createPrintAndSaveHistogram(wt_data, "ASCAT_smoothed_histograms/wt_hist_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "LOH histograms:"
lsl.createPrintAndSaveHistogram(wt_data, "ASCAT_smoothed_histograms/loh_hist_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "All Gain histograms:"
lsl.createPrintAndSaveHistogram(gain_data, "ASCAT_smoothed_histograms/gain_hist_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "3 Gain histograms:"
lsl.createPrintAndSaveHistogram(gain3_data, "ASCAT_smoothed_histograms/gain_hist3_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "4 Gain histograms:"
lsl.createPrintAndSaveHistogram(gain4_data, "ASCAT_smoothed_histograms/gain_hist4_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "5 Gain histograms:"
lsl.createPrintAndSaveHistogram(gain5_data, "ASCAT_smoothed_histograms/gain_hist5_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "6 Gain histograms:"
lsl.createPrintAndSaveHistogram(gain6_data, "ASCAT_smoothed_histograms/gain_hist6_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "Overall histogram:"
lsl.createPrintAndSaveHistogram(all_data, "ASCAT_smoothed_histograms/all_hist_a" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "00 histogram:"
lsl.createPrintAndSaveHistogram(data_00, "ASCAT_smoothed_histograms/00_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "01 histogram:"
lsl.createPrintAndSaveHistogram(data_01, "ASCAT_smoothed_histograms/01_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "02 histogram:"
lsl.createPrintAndSaveHistogram(data_02, "ASCAT_smoothed_histograms/02_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "03 histogram:"
lsl.createPrintAndSaveHistogram(data_03, "ASCAT_smoothed_histograms/03_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "11 histogram:"
lsl.createPrintAndSaveHistogram(data_11, "ASCAT_smoothed_histograms/11_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "12 histogram:"
lsl.createPrintAndSaveHistogram(data_12, "ASCAT_smoothed_histograms/12_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "13 histogram:"
lsl.createPrintAndSaveHistogram(data_13, "ASCAT_smoothed_histograms/13_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "22 histogram:"
lsl.createPrintAndSaveHistogram(data_22, "ASCAT_smoothed_histograms/22_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "23 histogram:"
lsl.createPrintAndSaveHistogram(data_23, "ASCAT_smoothed_histograms/23_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)

print "33 histogram:"
lsl.createPrintAndSaveHistogram(data_33, "ASCAT_smoothed_histograms/33_hist_" + rangestr + "all.txt", g_binwidth, axis=thisaxis)



print "Total number of points for range", rangestr, ": ", len(all_data)
