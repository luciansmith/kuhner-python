#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:50:09 2017

@author: lpsmith
"""

from __future__ import division
from os import walk

import numpy

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

flow_file = "FlowDataForChallenge.txt"

tag = "v6"

root_dir = "summary_stats/" + tag + "/"
ploidy_key = "ploidy"
purity_key = "purity"
psi_key = "psi"
goodness_key = "goodness"
unconstrained_key = "unconstrained"
diploid_key = "diploid"
tetraploid_key = "tetraploid"

ASCATvals = [ploidy_key, purity_key, psi_key, goodness_key]
constraints = [unconstrained_key, diploid_key, tetraploid_key]
outfile = "joint_summary_" + tag + ".txt"

files = {}

for ASCATval in ASCATvals:
    files[ASCATval] = {}
    for constraint in constraints:
        for (_, _, f) in walk(root_dir + ASCATval + "/" + constraint + "/"):
            files[ASCATval][constraint] = f

values = {}
patient_samples = set()

for ASCATval in files:
    values[ASCATval] = {}
    for constraint in files[ASCATval]:
        values[ASCATval][constraint] = {}
        for f in files[ASCATval][constraint]:
            infile = open(root_dir + ASCATval + "/" + constraint + "/" + f, "r")
            for line in infile:
                if line.find("x") != -1:
                    continue
                (id, val) = line.split()
                (patient, sample) = id.split("_")
                val = float(val)
                values[ASCATval][constraint][(patient, sample)] = val
                patient_samples.add((patient, sample))

patient_samples = sorted(patient_samples)

flow_data = {}
flowfile = open(flow_file, "r")
for line in flowfile:
    if line.find("Random") != -1:
        continue
    (patient, __, __, __, __, __, bxg2, a1Ploidy, percentA1, a2Ploidy, percentA2, __) = line.split("\t")
    if patient not in flow_data:
        flow_data[patient] = [0, 0, [], set()]
    if a1Ploidy == "":
        flow_data[patient][0] += 1
        flow_data[patient][2].append(0)
        flow_data[patient][3].add(2.0)
    else:
        flow_data[patient][3].add(float(a1Ploidy))
    if percentA1 != "":
        flow_data[patient][1] += 1
        percentA1 = float(percentA1)
        if (percentA2 != ""):
            percentA1 += float(percentA2)
        flow_data[patient][2].append(percentA1)
    if a2Ploidy != "":
        flow_data[patient][3].add(float(a2Ploidy))

for patient in flow_data:
    flow_data[patient][2] = numpy.average(flow_data[patient][2])

def hasCloseEntryIn(key, values):
    for value in values:
        if abs(key-value) < 0.2:
            return True
    return False


uncon_key = "unconstrained matches"
four_key = "exactly four"
#bpurity_key = "best purity"
dclose_key = "has close to diploid entry"
tclose_key = "has close to tetraploid entry"
ndreads = "diploid:non-diploid flows"
goodEv_key = "goodness diff"

evidence_types = [uncon_key, four_key, dclose_key, tclose_key, ndreads, goodEv_key]

evidence = {}
for etype in evidence_types:
    evidence[etype] = {}

for p_s in patient_samples:
    #Evidence that the unconstrained run matches one of the other runs
    if p_s in values[ploidy_key][unconstrained_key] and p_s in values[ploidy_key][diploid_key]:
        if values[ploidy_key][unconstrained_key][p_s] == values[ploidy_key][diploid_key][p_s]:
            evidence[uncon_key][p_s] = diploid_key
    if p_s in values[ploidy_key][unconstrained_key] and p_s in values[ploidy_key][tetraploid_key]:
        if values[ploidy_key][unconstrained_key][p_s] == values[ploidy_key][tetraploid_key][p_s]:
            evidence[uncon_key][p_s] = tetraploid_key

    #Evidence that the tetraploid result is exactly 4:
    if p_s in values[ploidy_key][tetraploid_key]:
        val = values[ploidy_key][tetraploid_key][p_s]
        if val < 4.1 and val > 3.9:
            evidence[four_key][p_s] = diploid_key
        else:
            evidence[four_key][p_s] = tetraploid_key

    #Evidence for the best purity:
#    uval = values[purity_key][unconstrained_key][p_s]
#    dval = 0
#    tval = 0
#    best = unconstrained_key
#    if p_s in values[purity_key][diploid_key]:
#        dval = values[purity_key][diploid_key][p_s]
#    if p_s in values[purity_key][tetraploid_key]:
#        tval = values[purity_key][tetraploid_key][p_s]
#    if dval >= uval and dval > tval:
#        best = diploid_key
#    elif tval >= uval and tval > dval:
#        best = tetraploid_key
#    elif tval >= uval and tval == dval:
#        best = "equal"
#    evidence[bpurity_key][p_s] = best

    #Evidence for particular ploidy showing up in flow data:
    if p_s[0] in flow_data:
        flow = flow_data[p_s[0]]
        if p_s in values[ploidy_key][diploid_key]:
            if hasCloseEntryIn(values[ploidy_key][diploid_key][p_s], flow[3]):
                evidence[dclose_key][p_s] = "True"
            else:
                evidence[dclose_key][p_s] = "False"
        if p_s in values[ploidy_key][tetraploid_key]:
            if hasCloseEntryIn(values[ploidy_key][tetraploid_key][p_s], flow[3]):
                evidence[tclose_key][p_s] = "True"
            else:
                evidence[tclose_key][p_s] = "False"
        evidence[ndreads][p_s] = str(flow[0]) + "::"+str(flow[1])
    else:
        evidence[ndreads][p_s] = "0::0"

    #Evidence: goodness of fit difference
    if p_s in values[goodness_key][tetraploid_key]:
        tetval = values[goodness_key][tetraploid_key][p_s]
        if p_s in values[goodness_key][diploid_key]:
            dipval = values[goodness_key][diploid_key][p_s]
            evidence[goodEv_key][p_s] = str(dipval - tetval)
        else:
            evidence[goodEv_key][p_s] = "NA"

out = open(outfile, "w")
out.write("patient\tsample")
for constraint in constraints:
    for ASCATval in ASCATvals:
        out.write("\t" + constraint + "_" + ASCATval)
for etype in evidence_types:
    out.write("\t" + etype)
out.write("\n")

for p_s in patient_samples:
    out.write(p_s[0] + "\t" + p_s[1])
    for constraint in constraints:
        for ASCATval in ASCATvals:
            if p_s in values[ASCATval][constraint]:
                out.write("\t" + str(values[ASCATval][constraint][p_s]))
            else:
                out.write("\tNA")
    for etype in evidence_types:
        if p_s in evidence[etype]:
            out.write("\t" + evidence[etype][p_s])
        else:
            out.write("\tNA")
    out.write("\n")
out.close()
