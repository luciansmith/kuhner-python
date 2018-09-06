#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:54:57 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

evidence_file = "calling_evidence.tsv"
outfile= "evidence_counts.tsv"

onlysomepatients = False
somepatients = ["572"]

def readEvidence():
    evidence = {}
    fcf = open(evidence_file, "r")
    for line in fcf:
        if "Patient" in line:
            continue
        (patient, sample, tomcall, VAFcat, flowratio, close_dip, close_tet, better_acc, acc_diff, NYGC_closer, goodness, Xcloser, finalcall) = line.rstrip().split("\t")
        if patient not in evidence:
            evidence[patient] = {}
        evidence[patient][sample] = {}
        evidence[patient][sample]["Tom's Partek Call"] = tomcall
        evidence[patient][sample]["2N VAF histogram category"] = VAFcat
        evidence[patient][sample]["Close diploid flow?"] = close_dip
        evidence[patient][sample]["Close tetraploid flow?"] = close_tet
        evidence[patient][sample]["Better accuracy"] = better_acc
        evidence[patient][sample]["NYGC closer to:"] = NYGC_closer
        evidence[patient][sample]["Goodness diff"] = goodness
        evidence[patient][sample]["Xiaohong closer"] = Xcloser
        evidence[patient][sample]["final"] = finalcall
    return evidence

def countCategories(evidence):
    counts = {}
    for patient in evidence:
        for sample in evidence[patient]:
            finalcall = evidence[patient][sample]["final"]
            if finalcall not in ("Diploid", "Tetraploid"):
                continue
            for category in evidence[patient][sample]:
                if category == "final":
                    continue
                if category not in counts:
                    counts[category] = {}
                call = evidence[patient][sample][category]
                if call not in counts[category]:
                    counts[category][call] = {}
                    counts[category][call]["Diploid"] = 0
                    counts[category][call]["Tetraploid"] = 0
                counts[category][call][finalcall] += 1
    return counts

def writeCounts(counts):
    outcounts = open(outfile, "w")
    outcounts.write("Evidence Category")
    outcounts.write("\tEvidence call")
    outcounts.write("\tnum_dip")
    outcounts.write("\tdip_running_total")
    outcounts.write("\tnum_tet")
    outcounts.write("\ttet_running_total")
    outcounts.write("\n")
    for category in counts:
        dip_total = 0
        tet_total = 0
        for call in counts[category]:
            outcounts.write(category)
            outcounts.write("\t" + call)
            outcounts.write("\t" + str(counts[category][call]["Diploid"]))
            dip_total += counts[category][call]["Diploid"]
            outcounts.write("\t" + str(dip_total))
            outcounts.write("\t" + str(counts[category][call]["Tetraploid"]))
            tet_total += counts[category][call]["Tetraploid"]
            outcounts.write("\t" + str(tet_total))
            outcounts.write("\n")
            

evidence = readEvidence()
counts = countCategories(evidence)
writeCounts(counts)