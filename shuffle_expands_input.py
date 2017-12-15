#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:22:45 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile

import lucianSNPLibrary as lsl
import random

expands_directory = "expands_full_input_alldupes_noLOH/"
shuffled_directory_core = "expands_full_input_shuffled_noLOH"

flist = []
filesets = {}
for (_, _, f) in walk(expands_directory):
    flist += f
    
seg_boundaries = (5796, 11592, 17388, 23184, 28981)
seg_boundaries = (55867, 643233, 7661420, 68082428)
bafs_by_seglengths = []
CNs_by_seglengths = []
for n in range(5):
    bafs_by_seglengths.append([])
    CNs_by_seglengths.append([])

def shuffleInParallel(list1, list2):
    if (len(list1) != len(list2)):
        print "Error: two lists are not the same length; cannot shuffle them in parallel."
        return
    for n in range(len(list1)-1):
        s = random.randint(n, len(list1)-1)
        t1 = list1[n]
        t2 = list2[n]
        list1[n] = list1[s]
        list2[n] = list2[s]
        list1[s] = t1
        list2[s] = t2
    
for f in flist:
    if f.find("BAF") == -1:
        continue
    baffile = open(expands_directory + f, "r")
    for line in baffile:
        if (line.find("chr") != -1):
            continue
        (chr, start, end, baf, PN_B) = line.split()
        seglength = int(end) - int(start)
        if (seglength <= seg_boundaries[0]):
            bafs_by_seglengths[0].append(baf)
        elif(seglength <= seg_boundaries[1]):
            bafs_by_seglengths[1].append(baf)
        elif(seglength <= seg_boundaries[2]):
            bafs_by_seglengths[2].append(baf)
        elif(seglength <= seg_boundaries[3]):
            bafs_by_seglengths[3].append(baf)
        else:
            bafs_by_seglengths[4].append(baf)
    
    CNname = f.replace("BAF", "CN")
    CNfile = open(expands_directory + CNname, "r")
    for line in CNfile:
        if (line.find("chr") != -1):
            continue
        (chr, start, end, CN) = line.split()
        seglength = int(end) - int(start)
        if (seglength <= seg_boundaries[0]):
            CNs_by_seglengths[0].append(CN)
        elif(seglength <= seg_boundaries[1]):
            CNs_by_seglengths[1].append(CN)
        elif(seglength <= seg_boundaries[2]):
            CNs_by_seglengths[2].append(CN)
        elif(seglength <= seg_boundaries[3]):
            CNs_by_seglengths[3].append(CN)
        else:
            CNs_by_seglengths[4].append(CN)
    if len(bafs_by_seglengths[0]) != len(CNs_by_seglengths[0]):
        print "Problem at file", f, "and", CNname
        random.nothing()

for sh in range(100):
    shuffled_directory = shuffled_directory_core + "_" + str(sh) + "/"
    for n in range(5):
        shuffleInParallel(bafs_by_seglengths[n], CNs_by_seglengths[n])
    print "Writing output for", shuffled_directory
    baf_indexes = [0, 0, 0, 0, 0]
    CN_indexes = [0, 0, 0, 0, 0]
    for f in flist:
        if f.find("BAF") == -1:
            continue
        #print "Processing file", f
        baffile = open(expands_directory + f, "r")
        bafout = open(shuffled_directory + f, "w")
        bafout.write("chr\tstartpos\tendpos\tAF_Tumor\tPN_B\n")
        for line in baffile:
            if (line.find("chr") != -1):
                continue
            (chr, start, end, baf, PN_B) = line.split()
            seglength = int(end) - int(start)
            if (seglength <= seg_boundaries[0]):
                baf = bafs_by_seglengths[0][baf_indexes[0]]
                baf_indexes[0] += 1
            elif (seglength <= seg_boundaries[1]):
                baf = bafs_by_seglengths[1][baf_indexes[1]]
                baf_indexes[1] += 1
            elif (seglength <= seg_boundaries[2]):
                baf = bafs_by_seglengths[2][baf_indexes[2]]
                baf_indexes[2] += 1
            elif (seglength <= seg_boundaries[3]):
                baf = bafs_by_seglengths[3][baf_indexes[3]]
                baf_indexes[3] += 1
            else:
                baf = bafs_by_seglengths[0][baf_indexes[4]]
                baf_indexes[4] += 1
            bafout.write(chr + "\t" + start + "\t" + end + "\t" + baf + "\t" + PN_B + "\n")
        bafout.close()

        CNname = f.replace("BAF", "CN")
        CNfile = open(expands_directory + CNname, "r")
        CNout = open(shuffled_directory + CNname, "w")
        CNout.write("chr\tstartpos\tendpos\tCN_Estimate\n")
        for line in CNfile:
            if (line.find("chr") != -1):
                continue
            (chr, start, end, CN) = line.split()
            seglength = int(end) - int(start)
            seglength = int(end) - int(start)
            if (seglength <= seg_boundaries[0]):
                CN = CNs_by_seglengths[0][CN_indexes[0]]
                CN_indexes[0] += 1
            elif (seglength <= seg_boundaries[1]):
                CN = CNs_by_seglengths[1][CN_indexes[1]]
                CN_indexes[1] += 1
            elif (seglength <= seg_boundaries[2]):
                CN = CNs_by_seglengths[2][CN_indexes[2]]
                CN_indexes[2] += 1
            elif (seglength <= seg_boundaries[3]):
                CN = CNs_by_seglengths[3][CN_indexes[3]]
                CN_indexes[3] += 1
            else:
                CN = CNs_by_seglengths[4][CN_indexes[4]]
                CN_indexes[4] += 1
            CNout.write(chr + "\t" + start + "\t" + end + "\t" + CN + "\n")
        CNout.close()
    
    