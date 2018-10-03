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
from os import mkdir
from os.path import isfile
from random import shuffle

import lucianSNPLibrary as lsl

gdir = "gamma_test_output/pASCAT_input_"
outfile = "purities_and_ploidies.tsv"

def getOddsAndPloidies():
    oddsfile = "calling_evidence_challenge_inc_odds.tsv"
    oddsAndPloidies = {}
    for line in open(oddsfile, "r"):
        if "Patient" in line:
            continue
        lvec = line.rstrip().split("\t")
        (patient, sample) = lvec[0:2]
        if patient not in oddsAndPloidies:
            oddsAndPloidies[patient] = {}
        #first check if the hutch made a decision about this sample:
        odds = float(lvec[-2])
        ploidy = lvec[-1]
        if ploidy == "Unknown":
            if odds >= 0.5:
                ploidy = "Diploid"
            else:
                ploidy = "Tetraploid"
        oddsAndPloidies[patient][sample] = (odds, ploidy)
    return oddsAndPloidies

def getPloidies(gdir, patient, ploidies):
    fname = gdir + patient + "_fcn_ascat_ploidy.txt"
    if not isfile(fname):
        return
    pfile = open(fname, "r")
    for line in pfile:
        if "x" in line:
            continue
        (__, sample, ploidy) = line.split('"')
        ploidy = float(ploidy)
        sample = sample.split("_")[1]
        ploidies[sample] = ploidy

def getPurities(gdir, patient, purities):
    fname = gdir + patient + "_fcn_ascat_cont.txt"
    if not isfile(fname):
        return
    pfile = open(gdir + patient + "_fcn_ascat_cont.txt", "r")
    for line in pfile:
        if "x" in line:
            continue
        (__, sample, purity) = line.split('"')
        purity = float(purity)
        sample = sample.split("_")[1]
        purities[sample] = purity

def writeHeaders(out):
    out.write("Patient")
    out.write("\tSample")
    out.write("\tOdds diploid better")
    out.write("\tPurity")
    out.write("\tPloidy")
    out.write("\n")

out = open(outfile, "w")
writeHeaders(out)

ploidies = {}
purities = {}
oddsAndPloidies = getOddsAndPloidies()
for patient in oddsAndPloidies:
    ploidies[patient] = {}
    purities[patient] = {}
    for constraint in ["diploid", "tetraploid"]:
        ploidies[patient][constraint] = {}
        purities[patient][constraint] = {}
    gamma = "g" + lsl.getGammaFor(patient)
    ddir = gdir + gamma + "/diploid/"
    tdir = gdir + gamma + "/tetraploid/"
    getPloidies(ddir, patient, ploidies[patient]["diploid"])
    getPurities(ddir, patient, purities[patient]["diploid"])
    getPloidies(tdir, patient, ploidies[patient]["tetraploid"])
    getPurities(tdir, patient, purities[patient]["tetraploid"])
    for sample in oddsAndPloidies[patient]:
        out.write(patient)
        out.write("\t" + sample)
        out.write("\t" + str(oddsAndPloidies[patient][sample][0]))
        if oddsAndPloidies[patient][sample][1] == "Diploid":
            out.write("\t" + str(purities[patient]["diploid"][sample]))
            out.write("\t" + str(ploidies[patient]["diploid"][sample]))
        elif oddsAndPloidies[patient][sample][1] == "Tetraploid":
            out.write("\t" + str(purities[patient]["tetraploid"][sample]))
            out.write("\t" + str(ploidies[patient]["tetraploid"][sample]))
        else:
            print("Unknown solution:", oddsAndPloidies[patient][sample][1], "for", patient, sample)
            assert(False)
        out.write("\n")

