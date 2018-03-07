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

analysis_dir = "gamma_test_output/analysis_compare/"
outdir = "best_analyses/"
gamma_list = ["100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000"]

onlysomepatients = False
somepatients = ["43"]

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

files = []
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
        (patient, sample, gamma, ploidy, TPn, FPn, UPn, TNn, FNn, UNn, NCn, Sn, TPl, FPl, UPl, TNl, FNl, UNl, NCl, Sl, __, __, __, __, __, __) = line.split()
        if TPn == "0" and FPn == "0" and UPn == "0":
            #This just failed?  I guess?
            continue
        if sample not in analysis_summaries:
            analysis_summaries[sample] = {}
        if ploidy not in analysis_summaries[sample]:
            analysis_summaries[sample][ploidy] = {}
        analysis_summaries[sample][ploidy][gamma] = {}
        analysis_summaries[sample][ploidy][gamma]["by_segment"] = (int(TPn) + int(TNn)) / (int(TPn) + int(FPn) + int(TNn) + int(FNn))
        analysis_summaries[sample][ploidy][gamma]["by_length"] =  (int(TPl) + int(TNl)) / (int(TPl) + int(FPl) + int(TNl) + int(FNl))
        if TPn == "0" and FPn == "0" and UPn == "0":
            #This means that ASCAT actually failed for this sample entirely, not that the 'true negatives' were wonderful.
            analysis_summaries[sample][ploidy][gamma]["by_segment"] = 0
            analysis_summaries[sample][ploidy][gamma]["by_length"] =  0
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
                best_out.write("\t" + best[sample][ploidy][segorlen][0])
                for gamma in gamma_list:
                    if ploidy not in analysis_summaries[sample] or gamma not in analysis_summaries[sample][ploidy]:
                        continue
                    if bestv[sample][ploidy][segorlen]==0 or analysis_summaries[sample][ploidy][gamma][segorlen]/bestv[sample]["overall"][segorlen] > 0.95:
                        best_out.write("\t" + gamma)
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
                for ploidy in ["diploid", "tetraploid"]:
                    for gamma in ["3000", "1000", "400", "250", "100"]:
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

#for segorlen in gamma_levels:
#    print("Histograms for", segorlen)
#    for gamma in gamma_levels[segorlen]:
#        print("Histogram for a gamma of", gamma)
#        lsl.createPrintAndSaveHistogram(gamma_levels[segorlen][gamma], "", 0.001)