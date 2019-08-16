#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 14:29:41 2019

@author: lpsmith
"""

tomfilename = "tom_final_DNA_calls.tsv"
tomcounts = "tom_counts.tsv"

counts = {}
calls = set()
for line in open(tomfilename, "r"):
    if "Call" in line:
        continue
    (call, patient, sample, nygc) = line.rstrip().split('\t')
    if patient not in counts:
        counts[patient] = {}
    if call not in counts[patient]:
        counts[patient][call] = 0
    counts[patient][call] += 1
    calls.add(call)

calls = list(calls)
calls.sort()

tc = open(tomcounts, "w")
tc.write("Patient")
for call in calls:
    tc.write("\t" + call)
tc.write("\n")

for patient in counts:
    tc.write(patient)
    for call in calls:
        if call in counts[patient]:
            tc.write("\t" + str(counts[patient][call]))
        else:
            tc.write("\t0")
    tc.write("\n")
tc.close()