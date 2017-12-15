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

print "Loading 2.5M and 1M SNPs."
m1labels, m1rev_labels = lsl.getSNPLabels(True, True)
m25labels, m25rev_labels = lsl.getSNPLabels(False, True)

all_out = open("probe_set_1_25_all.txt", "w")
overlap_out = open("probe_set_1_25_strict_overlap.txt", "w")
averaged_out = open("probe_set_1_25_averaged_overlap.txt", "w")
only_averaged_out = open("probe_set_1_25_only_averaged_overlap.txt", "w")

all_out.write("Name\tChr\tPosition\n")
overlap_out.write("Name\tChr\tPosition\n")
averaged_out.write("Name\tChr\tPosition\n")
only_averaged_out.write("Name\tChr\tPosition\n")

#The minimum distance between two SNPs allowed before averaging them:
min_avg_dist = 10000

m25_sorted = {}
m1_moved = {}
m25_used = {}

for label in m25labels:
    (chr, pos) = m25labels[label]
    if chr not in m25_sorted:
        m25_sorted[chr] = []
    m25_sorted[chr].append(int(pos))

doubled = 0
print "Writing out overlaps."
for label in m1labels:
    m1l = m1labels[label]
    if label in m25labels:
        (chr, pos) = m25labels[label]
        all_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
        overlap_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
        averaged_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
        m25_used[label] = (chr, pos)
#        if int(pos) in m25_sorted[chr]:
#            m25_sorted[chr].remove(int(pos))
    else:
        (chr, pos) = m1l
        all_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
        if m1l in m25rev_labels:
            m25label = m25rev_labels[m1l]
            overlap_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
            overlap_out.write(m25label + "\t" + chr + "\t" + str(pos) + "\n")
            averaged_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
            averaged_out.write(m25label + "\t" + chr + "\t" + str(pos) + "\n")
            m25_used[m25label] = (chr, pos)
#            if int(pos) in m25_sorted[chr]:
#                m25_sorted[chr].remove(int(pos))
        else:
            m1_moved[label] = m1l

overlap_out.close()

print "Writing 2.5M SNPs that weren't in 1M list."
for label in m25labels:
    if label in m1labels:
        continue
    (chr, pos) = m25labels[label]
    all_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")

all_out.close()

for chr in m25_sorted:
    print "Sorting chromosome", chr
    m25_sorted[chr].sort()

print "Averaging positions of close-enough SNPs from 1M and 2.5M"
for label in m1_moved:
    m1l = m1_moved[label]
    (chr, pos) = m1l
    if chr=="0":
        continue
    if pos=="0":
        continue
    pos = int(pos)
    candidate = -1000000000000
    for m25pos in m25_sorted[chr]:
        if abs(m25pos - pos) < abs(m25pos - candidate) and abs(m25pos - pos) < min_avg_dist:
            candidate = m25pos
        if m25pos - pos > 0:
            break
    m25l = (chr, str(candidate))
    if m25l in m25rev_labels:
        m25label = m25rev_labels[m25l]
        if m25label not in m25_used:
            avg_pos = int(math.floor((pos + candidate)/2))
            averaged_out.write(label + "\t" + chr + "\t" + str(avg_pos) + "\n")
            averaged_out.write(m25label + "\t" + chr + "\t" + str(avg_pos) + "\n")
            only_averaged_out.write(label + "\t" + chr + "\t" + str(avg_pos) + "\n")
            only_averaged_out.write(m25label + "\t" + chr + "\t" + str(avg_pos) + "\n")
            m25_used[m25label] = (chr, avg_pos)
        else:
            #We just add another option to the list
            (chr, pos) = m25_used[m25label]
            averaged_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
            only_averaged_out.write(label + "\t" + chr + "\t" + str(pos) + "\n")
#        m25_sorted[chr].remove(candidate)
#    else:
#        print "No candidate found for", label, m1l
        
             
averaged_out.close()
only_averaged_out.close()
        
#print "new m1's", len(new_m1)
#print "new m25's", len(new_m25)
for chr in m25_sorted:
    print "m25 leftovers for chromosome", chr, ":", len(m25_sorted[chr])
print "Next-best-SNPs found:", len(m1_moved)
print "Doubled originally, and now have separate SNPs:", doubled
        
