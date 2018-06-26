#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:50:41 2018

@author: lpsmith
"""
from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

analysis_dir = "analysis_compare/"
pASCAT_root = "gamma_test_output/pASCAT_input_g"
outdir = "best_analyses/"
#gamma_list = ["100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500"]#, "3000"]
gamma_list = ["500",]

onlysomepatients = True
somepatients = ["572"]

if not(path.isdir(outdir)):
    mkdir(outdir)

gamma_levels = {}
gamma_levels["by_segment"] = {}
gamma_levels["by_length"] = {}
for gamma in gamma_list:
    gamma_levels["by_segment"][gamma] = []
    gamma_levels["by_length"][gamma] = []

catch_out = open(outdir + "catch_close.tsv", "w")
catch_out.write("Patient")
catch_out.write("\tSample")
catch_out.write("\tSegOrLen")
catch_out.write("\tBest match")
catch_out.write("\tOrig best gamma")
catch_out.write("\tOrig best ploidy")
catch_out.write("\n")

def noResultsFor(patient, sample, gamma, ploidy):
    if ploidy=="Xiaohong":
        return False
    if ploidy=="overall":
        return False
#    file = pASCAT_root + gamma + "/" + ploidy + "/" + patient + "_" + sample + "_raw_segments.txt"
#    print(file)
    return not isfile(pASCAT_root + gamma + "/" + ploidy + "/" + patient + "_" + sample + "_raw_segments.txt")

def getBestPatientGamma(analysis_summaries):
    gamma_sums = {}
    for gamma in gamma_list:
        gamma_sums[gamma] = 0
        for sample in analysis_summaries:
            better_gamma = "None"
            if 'diploid' in analysis_summaries[sample] and gamma in analysis_summaries[sample]['diploid']:
                better_gamma = analysis_summaries[sample]['diploid'][gamma]['by_length']
            if 'tetraploid' in analysis_summaries[sample] and gamma in analysis_summaries[sample]['tetraploid']:
                t_gamma = analysis_summaries[sample]['tetraploid'][gamma]['by_length']
                if better_gamma=="None" or t_gamma > better_gamma:
                    better_gamma = t_gamma
            if better_gamma == "None" or gamma_sums[gamma] == "None":
                gamma_sums[gamma] = "None"
            else:
                gamma_sums[gamma] += better_gamma
    return gamma_sums

def findAndPrintBestGammas(all_gsums):
    global_bests = {}
    for patient in all_gsums:
        for gamma in all_gsums[patient]:
            if gamma not in global_bests:
                global_bests[gamma] = 0
            if global_bests[gamma] == "None" or all_gsums[patient][gamma] == "None":
                global_bests[gamma] = "None"
            else:
                global_bests[gamma] += all_gsums[patient][gamma]
    globals_out = open(outdir + "combined_overall.tsv", "w")
    globals_out.write("Patient")
    for gamma in gamma_list:
        globals_out.write("\tg" + gamma)
    globals_out.write("\n")
    globals_out.write("overall")
    for gamma in gamma_list:
        globals_out.write("\t" + str(global_bests[gamma]))
    globals_out.write("\n")
    for patient in all_gsums:
        globals_out.write(patient)
        for gamma in gamma_list:
            globals_out.write("\t" + str(all_gsums[patient][gamma]))
        globals_out.write("\n")

files = []

all_best = {}
all_close = {}

all_gsums = {}

for (__, __, f) in walk(analysis_dir):
    files += f
for f in files:
    if "analysis_overview" not in f:
        continue
    patient = f.split("_")[0]
    if (onlysomepatients and patient not in somepatients):
        continue
    analysis_summaries = {}
    
    analysis_file = open(analysis_dir + f, "r")
    for line in analysis_file:
        if "Patient" in line:
            continue
        (patient, sample, gamma, ploidy, TPn, FPn, UPn, TNn, FNn, UNn, NCn, Sn, TPl, FPl, UPl, TNl, FNl, UNl, NCl, Sl, n_acc, l_acc) = line.split()
#        if TPn == "0" and FPn == "0" and UPn == "0":
#            #This just failed?  I guess?
#            continue
#            print("No positives called for", patient, sample, gamma, ploidy)
        if noResultsFor(patient, sample, gamma, ploidy):
            print("Skipping", patient, sample, gamma, ploidy)
            continue
        if sample not in analysis_summaries:
            analysis_summaries[sample] = {}
        if ploidy not in analysis_summaries[sample]:
            analysis_summaries[sample][ploidy] = {}
        analysis_summaries[sample][ploidy][gamma] = {}
        if TPn=="0" and TNn=="0" and FPn=="0" and FNn=="0":
            analysis_summaries[sample][ploidy][gamma]["by_segment"] = "??"
        else:
            analysis_summaries[sample][ploidy][gamma]["by_segment"] = (int(TPn) + int(TNn)) / (int(TPn) + int(FPn) + int(TNn) + int(FNn))
        if TPl=="0" and TNl=="0" and FPl=="0" and FNl=="0":
            analysis_summaries[sample][ploidy][gamma]["by_length"] =  "??"
        else:
            analysis_summaries[sample][ploidy][gamma]["by_length"] =  (int(TPl) + int(TNl)) / (int(TPl) + int(FPl) +    int(TNl) + int(FNl))
#        if TPn == "0" and FPn == "0" and UPn == "0":
            #This means that ASCAT actually failed for this sample entirely, not that the 'true negatives' were wonderful.
#            analysis_summaries[sample][ploidy][gamma]["by_segment"] = 0
#            analysis_summaries[sample][ploidy][gamma]["by_length"] =  0
    analysis_file.close()
    best = {}
    bestv = {}
    for sample in analysis_summaries:
        best[sample] = {}
        bestv[sample] = {}
        best[sample]["overall"] = {}
        bestv[sample]["overall"] = {}
        for ploidy in analysis_summaries[sample]:
            best[sample][ploidy] = {}
            bestv[sample][ploidy] = {}
            for gamma in analysis_summaries[sample][ploidy]:
                for segorlen in analysis_summaries[sample][ploidy][gamma]:
                    if segorlen not in best[sample][ploidy]:
                        best[sample][ploidy][segorlen] = ("None", "None")
                        bestv[sample][ploidy][segorlen] = 0
                    if segorlen not in best[sample]["overall"]:
                        best[sample]["overall"][segorlen] = ("None", "None")
                        bestv[sample]["overall"][segorlen] = 0
                    if analysis_summaries[sample][ploidy][gamma][segorlen] > bestv[sample][ploidy][segorlen]:
                        bestv[sample][ploidy][segorlen] = analysis_summaries[sample][ploidy][gamma][segorlen]
                        best[sample][ploidy][segorlen] = (gamma, ploidy)
                    if ploidy != "Xiaohong":
                        if analysis_summaries[sample][ploidy][gamma][segorlen] > bestv[sample]["overall"][segorlen]:
                            bestv[sample]["overall"][segorlen] = analysis_summaries[sample][ploidy][gamma][segorlen]
                            best[sample]["overall"][segorlen] = (gamma, ploidy)
#                    if analysis_summaries[sample][ploidy][gamma][segorlen] < 0:
#                        print("Negative analysis level:", patient, sample, ploidy, segorlen, gamma)
    
    all_gsums[patient] = getBestPatientGamma(analysis_summaries)

    best_out = open(outdir + patient + "_best.tsv", "w")
    best_out.write("Patient")
    best_out.write("\tSample")
    best_out.write("\tConstraint")
    best_out.write("\tComparing")
    best_out.write("\tAccuracy")
    best_out.write("\tBest Gamma")
    best_out.write("\tClose Gammas")
    best_out.write("\n")
    for sample in best:
        for ploidy in best[sample]:
            for segorlen in best[sample][ploidy]:
                best_out.write(patient)
                best_out.write("\t" + sample)
                best_out.write("\t" + ploidy)
                best_out.write("\t" + segorlen)
                best_out.write("\t" + str(bestv[sample][ploidy][segorlen]))
                this_best = best[sample][ploidy][segorlen][0]
                if this_best not in all_best:
                    all_best[this_best] = 0
                all_best[this_best] += 1
                if this_best not in all_close:
                    all_close[this_best] = 0
                all_close[this_best] += 1
                best_out.write("\t" + this_best)
                for gamma in gamma_list:
                    if ploidy not in analysis_summaries[sample] or gamma not in analysis_summaries[sample][ploidy]:
                        continue
                    if bestv[sample][ploidy][segorlen]==0 or analysis_summaries[sample][ploidy][gamma][segorlen]/bestv[sample]["overall"][segorlen] > 0.95:
                        best_out.write("\t" + gamma)
                    if gamma not in all_close:
                        all_close[gamma] = 0
                    all_close[gamma] += 1
                    if bestv[sample][ploidy][segorlen] == 0:
                        gamma_levels[segorlen][gamma].append(1)
                    elif analysis_summaries[sample][ploidy][gamma][segorlen] < 0:
                        gamma_levels[segorlen][gamma].append(0)
                    else:
                        gamma_levels[segorlen][gamma].append(analysis_summaries[sample][ploidy][gamma][segorlen]/bestv[sample][ploidy][segorlen])
                best_out.write("\n")

    for sample in best:
        for segorlen in best[sample]["overall"]:
            closestmatch = 0
            matches = "No Match"
            if bestv[sample]["overall"][segorlen]==0:
                matches = "Best match was zero."
            else:
                for ploidy in ["diploid", "tetraploid", "eight"]:
                    for gamma in gamma_list:
                        if ploidy in analysis_summaries[sample] and gamma in analysis_summaries[sample][ploidy]:
                            match = analysis_summaries[sample][ploidy][gamma][segorlen]/bestv[sample]["overall"][segorlen]
                            if match > closestmatch:
                                closestmatch = match
                matches = "Best Match:\t" + str(closestmatch)
            catch_out.write(patient)
            catch_out.write("\t" + sample)
            catch_out.write("\t" + segorlen)
            catch_out.write("\t" + str(closestmatch))
            catch_out.write("\t" + best[sample]["overall"][segorlen][0])
            catch_out.write("\t" + best[sample]["overall"][segorlen][1])
            catch_out.write("\n")
    
    best_out.close()

findAndPrintBestGammas(all_gsums)

print("All best:", all_best)
print("All close:", all_close)
#for segorlen in gamma_levels:
#    print("Histograms for", segorlen)
#    for gamma in gamma_levels[segorlen]:
#        print("Histogram for a gamma of", gamma)
#        lsl.createPrintAndSaveHistogram(gamma_levels[segorlen][gamma], "", 0.001)