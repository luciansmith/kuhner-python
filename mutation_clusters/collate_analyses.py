#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 11:52:00 2019

Collect the results of 'classify.py' (the 'analysis' files).

@author: lpsmith
"""

import os
from os import path
from os import mkdir
import lucianSNPLibrary as lsl
import matplotlib.pyplot as plt

somepatientsonly = False
somepatients = ["1005"]

analysisdir = "analysis/"
outdir = "collate/"

if not(path.isdir(outdir)):
    mkdir(outdir)

def readProgressorOrNot():
    progressorMap = {}
    pfile = open("../Patient_status.txt", "r")
    for line in pfile:
        if "Patient" in line:
            continue
        (patient, pstat, sex) = line.rstrip().split()
        progressorMap[patient] = pstat
    return progressorMap


def readGenomeDoubled():
    doubledMap = {}
    pfile = open("../calling_evidence_challenge_inc_odds.tsv", "r")
    for line in pfile:
        if "Patient" in line:
            continue
        lvec = line.rstrip().split()
        sample = lvec[1]
        call = lvec[-1]
        if call=="Unknown":
            odds = float(lvec[-2])
            if odds<.5:
                call = "Tetraploid"
            else:
                call = "Diploid"
        doubledMap[sample] = call
    return doubledMap


afiles = []
for root, dirs, files in os.walk(analysisdir):
  for file in files:
    if file.endswith("_analysis.txt"):
      afiles.append(file)


progressorMap = readProgressorOrNot()
doubledMap = readGenomeDoubled()

summary = {}
for file in afiles:
    (patient, sample, A, B, __) = file.split("_")
    if (A, B) not in summary:
        summary[(A, B)] = {}
    if patient not in summary[(A, B)]:
        summary[(A, B)][patient] = {}
    if sample not in summary[(A, B)][patient]:
        summary[(A, B)][patient][sample] = {}
    afile = open(analysisdir + file, "r")
    for line in afile:
        lvec = line.rstrip().split()
        if lvec[0] == "total":
            summary[(A, B)][patient][sample]["total"] = int(lvec[1])
        if lvec[0] == "cis/trans":
            summary[(A, B)][patient][sample]["ct_ratio"] = float(lvec[2])
        if lvec[0] == "cis":
            summary[(A, B)][patient][sample]["cis"] = int(lvec[1])
        if lvec[0] == "nested":
            summary[(A, B)][patient][sample]["nested"] = int(lvec[1])
        if lvec[0] == "trans":
            summary[(A, B)][patient][sample]["trans"] = int(lvec[1])
        if lvec[0] == "missed":
            summary[(A, B)][patient][sample]["missed"] = int(lvec[1])
        if lvec[0] == "wt":
            summary[(A, B)][patient][sample]["wt"] = int(lvec[1])
        if lvec[0] == "fourgamete":
            summary[(A, B)][patient][sample]["fourgamete"] = int(lvec[1])
        if lvec[0] == "anomaly":
            summary[(A, B)][patient][sample]["anomaly"] = int(lvec[1])
        if lvec[0] == "noreads":
            summary[(A, B)][patient][sample]["noreads"] = int(lvec[1])
        if lvec[0] == "toosmall":
            summary[(A, B)][patient][sample]["toosmall"] = int(lvec[1])

#for (A, B) in summary:
#    print("Summarizing over all CNs of", A, ",", B, ".")
#    hist = []
#    for patient in summary[(A, B)]:
#        for sample in summary[(A, B)][patient]:
#            hist.append(summary[(A, B)][patient][sample]["ct_ratio"])
#    lsl.createPrintAndSaveHistogram(hist, A + "_" + B + "_ct_ratio", 0.1, xdata="Cis/Trans Ratio")

cises = {}
transes = {}
for p in ("all", "Prog", "NP", "Diploid", "Tetraploid", "2_2"):
    cises[p] = []
    transes[p] = []

for patient in summary['1', '1']:
    prog = progressorMap[patient]
    for sample in summary['1', '1'][patient]:
        doubled = doubledMap[sample]
        cises["all"].append(summary['1', '1'][patient][sample]["cis"])
        transes["all"].append(summary['1', '1'][patient][sample]["trans"])
#        try:
#            cises["2_2"].append(summary['2', '2'][patient][sample]["cis"])
#            transes["2_2"].append(summary['2', '2'][patient][sample]["trans"])
#        except:
#            pass
        cises[prog].append(summary['1', '1'][patient][sample]["cis"])
        transes[prog].append(summary['1', '1'][patient][sample]["trans"])
        if doubled=="Diploid":
            cises[doubled].append(summary['1', '1'][patient][sample]["cis"])
            transes[doubled].append(summary['1', '1'][patient][sample]["trans"])
#        if doubled=="Tetraploid":
#            cises[doubled].append(summary['2', '2'][patient][sample]["cis"])
#            transes[doubled].append(summary['2', '2'][patient][sample]["trans"])

for patient in summary['2', '2']:
    prog = progressorMap[patient]
    for sample in summary['2', '2'][patient]:
        doubled = doubledMap[sample]
        cises["2_2"].append(summary['2', '2'][patient][sample]["cis"])
        transes["2_2"].append(summary['2', '2'][patient][sample]["trans"])
        if doubled=="Tetraploid":
            cises[doubled].append(summary['2', '2'][patient][sample]["cis"])
            transes[doubled].append(summary['2', '2'][patient][sample]["trans"])


print("Progressors (blue) vs. non-progressors (orange), for all cis vs. trans for CN ratios of 1,1")
plt.xlabel('cis')
plt.ylabel('trans')
plt.scatter(cises["Prog"], transes["Prog"])
plt.scatter(cises["NP"], transes["NP"])
plt.plot((0,700),(0,700))
plt.show()
plt.close()

print("cis vs. trans:  Blue: Diploid samples with CN ratios of 1,1.  Orange: Tetraploid samples with CN ratios of 2,2")
plt.xlabel('cis')
plt.ylabel('trans')
plt.scatter(cises["Diploid"], transes["Diploid"])
plt.scatter(cises["Tetraploid"], transes["Tetraploid"])
plt.plot((0,700),(0,700))
plt.show()
plt.close()

print("All cis vs. trans for CN ratios of 1,1 (blue) or 2,2 (orange)")
plt.xlabel('cis')
plt.ylabel('trans')
plt.scatter(cises["all"], transes["all"])
plt.scatter(cises["2_2"], transes["2_2"])
plt.plot((0,700),(0,700))
plt.show()
plt.close()

