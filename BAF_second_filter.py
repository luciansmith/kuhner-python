q#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 15:06:47 2016

@author: lpsmith
"""

import lucianSNPLibrary as lsl
from os import walk
from os import mkdir
from os import path
import string

tolerance = (0.4, 0.65)

onepatientonly = False
onepatient = "997"

tag = "_25M_only"
BAF_input = "BAF_first_filtered_data" + tag + "/"
BAF_output = "BAF_filtered_data" + tag + "_" + str(int(tolerance[0]*100)) + "_" + str(int(tolerance[1]*100)) + "/"
if not(path.isdir(BAF_output)):
    mkdir(BAF_output)

#filenames = ["omni15_paired_copynumber_baselineAdjusted_37_RB2012_2015.txt"]
filenames = []
for (_, _, f) in walk(BAF_input):
    filenames += f
    break

wt_files = set()

for file in filenames:
    if "BLD" in file or "gastric" in file:
        wt_files.add(file)
        continue
#    if (string.find(file, "N") != -1):
#        wt_files.add(file)
#        continue


# Patty says that we should keep the 'BLD2' samples and *not* the 'BLD' samples, as the former were run to replaced bad versions of the latter.
oldBLDs = set()
for f in wt_files:
    bld2 = f.find("BLD2")
    if (bld2 != -1):
        bld = f.replace("BLD2", "BLD")
        oldBLDs.add(bld)

for oldbld in oldBLDs:
    wt_files.remove(oldbld)

wt_data = {}
for file in wt_files:
    filebits = file.rstrip().split("_")
    patient = filebits[0]
    if onepatientonly and patient != onepatient:
        continue
    BAFfile = open(BAF_input + file,"r")
    wt_data[patient] = set()
    print("Checking WT BAFs for patient", patient)
    for line in BAFfile:
        if "BAF" in line:
            continue
        (id, chr, pos, baf) = line.rstrip().split("\t")
        if "cnvi" in id:
            continue
        try:
            baf = float(baf)
        except:
            #Unknown BAF value; probably '?'
            continue
        if (baf > tolerance[0]) and (baf < tolerance[1]):
            wt_data[patient].add(id)

for file in filenames:
    if "BLD" in file or "gastric" in file:
        continue
    filebits = file.rstrip().split("_")
    if (len(filebits) != 3):
        continue
    patient = filebits[0]
    if onepatientonly and patient != onepatient:
        continue
    sample = filebits[1]
    if (patient not in wt_data):
        print("No WT data for patient", patient,": skipping.")
        continue
    print("Filtering BAF output for patient", patient, "sample", sample, ".")
    BAFfile = open(BAF_input + file,"r")
    outfile = open(BAF_output + file, "w")
    outfile.write("SNPid\tchr\tpos\tBAF\n")
    for line in BAFfile:
        if (line.find("BAF") != -1):
            continue
        (id, chr, pos, baf) = line.rstrip().split("\t")
        if (id in wt_data[patient]):
            outfile.write(line)
    outfile.close()

