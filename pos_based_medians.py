#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

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

BAF_dirs = {}
BAF_dirs["1M"] = "BAF_first_filtered_data_1M_only/"
BAF_dirs["25M"] = "BAF_first_filtered_data_25M_only/"


medians = {}
medians["1M"] = {}
medians["25M"] = {}


for arraytype in BAF_dirs:
    bafdir = BAF_dirs[arraytype]
    files = []
    for (__, __, f) in walk(bafdir):
        files += f
    for f in files:
        if "BAF.txt" not in f:
            continue
        mlist = medians[arraytype]
        for line in open(bafdir + f, "r"):
            if "SNPid" in line:
                continue
            (id, chr, pos, baf) = line.rstrip().split()
            try:
                baf = float(baf)
                pos = int(pos)
                chr = int(chr)
            except:
                continue
            if (baf<0.4 or baf>.65):
                continue
            if chr not in mlist:
                mlist[chr] = {}
            if pos not in mlist[chr]:
                mlist[chr][pos] = []
            mlist[chr][pos].append(baf)
    
    mfile = open("medians_" + arraytype + ".tsv", "w")
    mfile.write("Chr\tpos\tmedian BAF\n")
    
    mlistkeys = list(mlist.keys())
    mlistkeys.sort();
    for chr in mlistkeys:
        poskeys = list(mlist[chr].keys())
        poskeys.sort()
        for pos in poskeys:
            mfile.write(str(chr))
            mfile.write("\t" + str(pos))
            mfile.write("\t" + str(numpy.median(mlist[chr][pos])))
            mfile.write("\n")
    mfile.close()
