#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:11:30 2017

@author: lpsmith
"""

from __future__ import division
import lucianSNPLibrary as lsl
from os import walk
import math
#from os import mkdir
#from os import path
#import string

labels = {}
print("Loading 2.5M v1.3, 2.5M v1.2, and 1M SNPs.")
labels["25Mv13"], __ = lsl.getSNPLabels2_5M(True)
labels["25Mv12"], __ = lsl.getSNPLabels2_5Mv12(True)
labels["1M"], __ = lsl.getSNPLabels1M(True)

all_out = open("probe_set_1_25_all.txt", "w")
all_out.write("Name\tChr\tPosition\n")

used_labels = set()

for probeset in ["25Mv13", "25Mv12", "1M"]:
    print("Starting probeset", probeset)
    for label in labels[probeset]:
        if label in used_labels:
            continue
        used_labels.add(label)
        (chr, pos) = labels[probeset][label]
        all_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")

all_out.close()

