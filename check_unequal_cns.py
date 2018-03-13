#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:43:10 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

import lucianSNPLibrary as lsl

CNdir = "nonintegerCNs/"

onlysomepatients = False
somepatients = ["997"]

firstpatients = ["17", "42", "55", "59", "74", "43", "184", "163", "396", "1047"]

files = []
firstmore = 0
secondmore = 0
firstmore_fromeq = 0
secondmore_fromeq = 0
firstmore_zeroB = 0
zeroB = 0
zeroB_zeroA = 0
total = 0
for (__, __, f) in walk(CNdir):
    files += f
for f in files:
    if "nonint" not in f:
        continue
    (patient, sample, ploidy) = f.split("_")[0:3]
    #print(patient, sample, ploidy, ":")
    cnfile = open(CNdir + f, "r")
    for line in cnfile:
        if "patient" in line:
            continue
        (patient, sample, chr, start, end, rawA, rawB, intA, intB) = line.split()
        if intA == "NA" or intB == "NA":
            continue
        total += 1
        intA = int(intA)
        intB = int(intB)
        if rawB == "0":
            zeroB += 1
            if rawA == "0" or intA == 0:
                zeroB_zeroA += 1
        if rawA == rawB and intA != intB:
            if intA > intB:
                firstmore_fromeq += 1
            else:
                secondmore_fromeq += 1
                print(line)
        elif intA != intB:
            if intA > intB:
                firstmore += 1
                if rawB == "0":
                    firstmore_zeroB += 1
                else:
                    print(line)
            else:
                secondmore += 1

print("First more:", firstmore)
print("Second more:", secondmore)
print("First more from equal raw values:", firstmore_fromeq)
print("Second more from equal raw values:", secondmore_fromeq)
print("Total segments:", total)
