#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 13:49:04 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
import random
import re
from os import path
from os import mkdir

import lucianSNPLibrary as lsl

tag = "_g1000_unconstrained/"

tree_input = "joint_segmentation_nodupes" + tag
beast_output = "BEAST" + tag

if not(path.isdir(beast_output + "/")):
    mkdir(beast_output + "/")

segment_files = []
for (_, _, f) in walk(tree_input):
    segment_files += f

for f in segment_files:
    if f.find("characters") == -1:
        continue
    patient = f.split("_")[0]
    segfile = open(tree_input + f, "r")
    samples = segfile.readline().split()[5:]
    A_out = open(beast_output + patient + "_phased_allA.txt", "w")
    B_out = open(beast_output + patient + "_phased_allB.txt", "w")
    A_out.write("Chr\tStart\tEnd\tnLogR\tnBAF")
    B_out.write("Chr\tStart\tEnd\tnLogR\tnBAF")
    for sample in samples:
        A_out.write("\t" + sample)
        B_out.write("\t" + sample)
    A_out.write("\n")
    B_out.write("\n")
    assignments = {}
    for line in segfile:
        dataline = line.split()
        (chr, start, end, nbaf, nCN) = dataline[0:5]
        A_out.write(chr + "\t" + start + "\t" + end + "\t" + nCN + "\t" + nbaf)
        B_out.write(chr + "\t" + start + "\t" + end + "\t" + nCN + "\t" + nbaf)
        for sample_data in dataline[5:]:
            (label1, data1, __, label2, data2)  = sample_data.split(":")
            label1 = re.sub('[?]', '', label1)
            label2 = re.sub('[?]', '', label2)
            if (label1.find("S") != -1):
                tmp1 = label1
                tmp2 = data1
                label1 = label2
                data1 = data2
                label2 = tmp1
                data2 = tmp2
            if label1 not in assignments:
                assignments[label1] = random.choice([True, False])
            if (assignments[label1]):
                A_out.write("\t" + data1)
                B_out.write("\t" + data2)
            else:
                A_out.write("\t" + data2)
                B_out.write("\t" + data1)
        A_out.write("\n")
        B_out.write("\n")
    A_out.close()
    B_out.close()

