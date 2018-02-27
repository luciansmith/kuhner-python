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
gamma_list = ["0", "50", "100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000"]

onlysomepatients = False
somepatients = ["43"]

if not(path.isdir(outdir)):
    mkdir(outdir)

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
        if sample not in analysis_summaries:
            analysis_summaries[sample] = {}
        if ploidy not in analysis_summaries[sample]:
            analysis_summaries[sample][ploidy] = {}
        analysis_summaries[sample][ploidy][gamma] = {}
        analysis_summaries[sample][ploidy][gamma]["by_segment"] = int(TPn) - int(FPn) + int(TNn) - int(FNn)
        analysis_summaries[sample][ploidy][gamma]["by_length"] =  int(TPl) - int(FPl) + int(TNl) - int(FNl)
    analysis_file.close()
    best = {}
    bestv = {}
    for sample in analysis_summaries:
        best[sample] = {}
        bestv[sample] = {}
        for ploidy in analysis_summaries[sample]:
            best[sample][ploidy] = {}
            bestv[sample][ploidy] = {}
            for gamma in analysis_summaries[sample][ploidy]:
                for segorlen in analysis_summaries[sample][ploidy][gamma]:
                    if segorlen not in best[sample][ploidy]:
                        best[sample][ploidy][segorlen] = ""
                        bestv[sample][ploidy][segorlen] = 0
                    if analysis_summaries[sample][ploidy][gamma][segorlen] > bestv[sample][ploidy][segorlen]:
                        bestv[sample][ploidy][segorlen] = analysis_summaries[sample][ploidy][gamma][segorlen]
                        best[sample][ploidy][segorlen] = gamma


    best_out = open(outdir + patient + "_best.tsv", "w")
    best_out.write("Patient")
    best_out.write("\tSample")
    best_out.write("\tConstraint")
    best_out.write("\tComparing")
    best_out.write("\tGood-Bad")
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
                best_out.write("\t" + best[sample][ploidy][segorlen])
                for gamma in analysis_summaries[sample][ploidy]:
                    if bestv[sample][ploidy][segorlen]==0 or analysis_summaries[sample][ploidy][gamma][segorlen]/bestv[sample][ploidy][segorlen] > 0.99:
                        best_out.write("\t" + gamma)
                best_out.write("\n")
    
    
    best_out.close()
