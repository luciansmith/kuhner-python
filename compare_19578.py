#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 13:54:45 2018

@author: lpsmith
"""
from __future__ import division
import lucianSNPLibrary as lsl

cnfiles = ["CN_filtered_data_1M_only/391_19578_copynumber_all.txt", "CN_filtered_data_Pilot/391_19578_copynumber_all.txt"]
baffiles = ["BAF_first_filtered_data_1M_only/391_19578_BAF.txt", "BAF_first_filtered_data_Pilot/391_19578_BAF.txt"]

cndata = {}
for cnfile in cnfiles:
    for line in open(cnfile, "r"):
        if "SNPid" in line:
            continue
        (snpid, chr, pos, l2r) = line.split()
        if snpid not in cndata:
            cndata[snpid] = [l2r]
        else:
            cndata[snpid].append(l2r)

bafdata = {}
for baffile in baffiles:
    for line in open(baffile, "r"):
        if "SNPid" in line:
            continue
        (snpid, chr, pos, baf) = line.split()
        if snpid not in bafdata:
            bafdata[snpid] = [baf]
        else:
            bafdata[snpid].append(baf)

cndiffs = []
bafdiffs = []
cns = [[], []]
bafs = [[], []]
for snpid in cndata:
    if len(cndata[snpid]) == 2:
        try:
            l2r_a = float(cndata[snpid][0])
            l2r_b = float(cndata[snpid][1])
            cndiffs.append(abs(l2r_a - l2r_b+0.32260464413499457))
            cns[0].append(l2r_a)
            cns[1].append(l2r_b)
        except:
#            print("not float values:", cndata[snpid])
            continue

for snpid in bafdata:
    if len(bafdata[snpid]) == 2:
        try:
            baf_a = float(bafdata[snpid][0])
            baf_b = float(bafdata[snpid][1])
            bafdiffs.append(abs(baf_a - baf_b))
            bafs[0].append(baf_a)
            bafs[1].append(baf_b)
        except:
#            print("not float values:", bafdata[snpid])
            continue

#print("CN diffs:")
#lsl.createPrintAndSaveHistogram(cndiffs, "", 0.1)

#print("BAF diffs:")
#lsl.createPrintAndSaveHistogram(bafdiffs, "", 0.01)




import numpy

print("Median BAF difference:", numpy.median(bafdiffs))
print("Median CN difference:", numpy.median(cndiffs))
print("Average BAF difference:", numpy.average(bafdiffs))
print("Average CN difference:", numpy.average(cndiffs))
print("CN differences >1:", sum(1 for x in cndiffs if x>1))
print("Total CN differences:", len(cndiffs))
print(437/688081)
print("Total BAF differences > 0.2:", sum(1 for x in bafdiffs if x>0.2))
print("Total BAF differences:", len(bafdiffs))
print(1550/688529)
