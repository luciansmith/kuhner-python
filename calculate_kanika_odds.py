# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 09:20:30 2019

@author: Lucian
"""

from __future__ import division
from os import walk
import matplotlib.pyplot as plt


kanikadir = "snv_union_tables_per_patient/"
sumfile = "nalt_calling_odds.tsv"


kanika_files = []
for __, _, files in walk(kanikadir):
    kanika_files += files

calling_odds = {}

for file in kanika_files:
    if "snv.union" not in file:
        continue
    patient = file.split(".")[0]
    allpos = {}
    badpos = set()
    for line in open(kanikadir + file, "r"):
        lvec = line.split("\t")
        if "CHROM" in line:
            headers = lvec
            continue
        (chrom, pos, ref, alt) = lvec[0:4]
        maxcallers = int(lvec[-1])
        if chrom=="X":
            continue
        if maxcallers < 2:
            continue
        if (chrom, pos) in allpos:
            badpos.add((chrom, pos))
        allpos[(chrom, pos)] = lvec[8:len(lvec)-3]
    for chrpos in badpos:
        del allpos[chrpos]
    for chrpos in allpos:
        vec = allpos[chrpos]
        for i in range(0, len(vec), 3):
            (VAF, nreads, ncalls) = vec[i:i+3]
            VAF = float(VAF)
            nreads = int(nreads)
            ncalls = int(ncalls)
            nalt = round(VAF*nreads)
            if nalt not in calling_odds:
                calling_odds[nalt] = [0, 0]
            if ncalls > 1:
                calling_odds[nalt][0] += 1
            calling_odds[nalt][1] += 1

nalts = list(calling_odds.keys())
nalts.sort()

x = []
y = []
summary = open(sumfile, "w")
summary.write("nAlt\tnTwoPlus\tnTot\tOdds\n")
for nalt in nalts:
    x.append(nalt)
    y.append(calling_odds[nalt][0]/calling_odds[nalt][1])
    summary.write(str(nalt))
    summary.write("\t" + str(calling_odds[nalt][0]))
    summary.write("\t" + str(calling_odds[nalt][1]))
    summary.write("\t" + str(calling_odds[nalt][0]/calling_odds[nalt][1]))
    summary.write("\n")
    if nalt <= 20:
        print(str(nalt) + "\t" + str(calling_odds[nalt][0]) + "\t" + str(calling_odds[nalt][1]) + "\t" + str(calling_odds[nalt][0]/calling_odds[nalt][1]))
summary.close()

plt.plot(x,y)