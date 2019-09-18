# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 15:54:51 2019

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
import csv
import ete3

import lucianSNPLibrary as lsl

delfile = "check_for_deletions.tsv"
CNdir = "noninteger_processed_CNs/"
outfile = "found_deletions_or_not.txt"


dels = {}

onlysomepatients = False
somepatients = ["1005"]


for line in open(delfile, "r"):
    if "Patient" in line:
        continue
    (patient, sample, chrom, pos, gene, incong) = line.rstrip().split()
    if onlysomepatients and patient not in somepatients:
        continue
    if patient not in dels:
        dels[patient] = {}
    if sample not in dels[patient]:
        dels[patient][sample] = {}
    if chrom not in dels[patient][sample]:
        dels[patient][sample][chrom] = {}
    dels[patient][sample][chrom][(int(pos), gene, incong)] = False

for patient in dels:
    for sample in dels[patient]:
        ploidy = lsl.getBestPloidyFor(patient, sample)
        filename = patient + "_" + sample + "_g500_" + ploidy + "_nonint_CNs.txt"
        if not isfile(CNdir + filename):
            filename = patient + "_" + sample + "_g550_" + ploidy + "_nonint_CNs.txt"
        assert(isfile(CNdir + filename))
        print("checking file", filename)
        for line in open(CNdir + filename, "r"):
            if "patient" in line:
                continue
            (patient, sample, chrom, start, end, __, __, intA, intB) = line.rstrip().split()
            if chrom in dels[patient][sample]:
                if intA=="NA" or intB=="NA":
                    continue
                intA = int(intA)
                intB = int(intB)
                start = int(start)
                end = int(end)
                if intA==0 or intB==0:
                    for pgi in dels[patient][sample][chrom]:
                        if pgi[0] >= start and pgi[0] <= end:
                                dels[patient][sample][chrom][pgi] = True

delfound = open(outfile, "w")
delfound.write("Patient\tSample\tchrom\tpos\tgene\tincongruous?\tFound Deletion?\n")


for patient in dels:
    for sample in dels[patient]:
        for chrom in dels[patient][sample]:
            for pgi in dels[patient][sample][chrom]:
                delfound.write(patient)
                delfound.write("\t" + sample)
                delfound.write("\t" + chrom)
                delfound.write("\t" + str(pgi[0]))
                delfound.write("\t" + pgi[1])
                delfound.write("\t" + pgi[2])
                delfound.write("\t" + str(dels[patient][sample][chrom][pgi]))
                delfound.write("\n")
delfound.close()    