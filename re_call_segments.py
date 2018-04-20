#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 15:18:36 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
import shutil

import numpy
import lucianSNPLibrary as lsl



nonint_dir = "nonintegerCNs/"
balance_dir = "balanced_calls/"

outdir = "noninteger_processed_CNs/"
jamboreedir = "jamboree_files/"

onlysomepatients = False
somepatients = ["55"]

if not path.isdir(outdir):
    mkdir(outdir)

if not path.isdir(jamboreedir + outdir):
    mkdir(jamboreedir + outdir)


def getNonIntegerCalls(f):
    nonints = {}
    nonint = open(nonint_dir + f, "r")
    for line in nonint:
        if "patient" in line:
            continue
        (patient, sample, chr, start, end, rawA, rawB, intA, intB) = line.split()
        try:
            rawA = float(rawA)
            rawB = float(rawB)
        except:
            rawA = numpy.nan
            rawB = numpy.nan
        if chr not in nonints:
            nonints[chr] = []
        nonints[chr].append([start, end, rawA, rawB, intA, intB])
    return nonints

def modifyDifficultCalls(nonints, bal_calls):
    #For anything where rawA = rawB (ish), and intA != intB, we change it 
    # so that intA = intB, if the 'balanced' algorithm says it's supposed to match.
    for chr in nonints:
        for seg in nonints[chr]:
            (start, end, rawA, rawB, intA, intB) = seg
            start = int(start)
            end = int(end)
            if intA != intB and abs(rawA - rawB) < 0.1:
                for bseg in bal_calls[chr]:
                    if bseg[0] <= start and bseg[1] >= end:
                        if bal_calls[chr][bseg] == "Balanced":
                            closeint = round((rawA + rawB)/2)
                            seg[4] = str(closeint)
                            seg[5] = str(closeint)
                        #else, we leave it as it was

def writeNewNonInts(f, nonints):
    nonint = open(outdir + f, "w")
    (patient, sample, gamma, ploidy) = f.split("_")[0:4]
    nonint.write("patient")
    nonint.write("\tbiopsy")
    nonint.write("\tchrom")
    nonint.write("\tsegstart")
    nonint.write("\tsegend")
    nonint.write("\trawA")
    nonint.write("\trawB")
    nonint.write("\tintA")
    nonint.write("\tintB")
    nonint.write("\n")
    for chr in nonints:
        for seg in nonints[chr]:
            nonint.write(patient)
            nonint.write("\t" + sample)
            nonint.write("\t" + chr)
            nonint.write("\t" + seg[0])
            nonint.write("\t" + seg[1])
            if numpy.isnan(seg[2] or numpy.isnan(seg[3])):
                nonint.write("\tNA\tNA")
            else:
                nonint.write("\t" + str(seg[2]))
                nonint.write("\t" + str(seg[3]))
            nonint.write("\t" + seg[4])
            nonint.write("\t" + seg[5])
            nonint.write("\n")
    nonint.close()

def copyToJamboreeIfWGS(f):
    (patient, sample, gamma, ploidy) = f.split("_")[0:4]
    if int(sample) < 23341 and sample != "19578":
        return
    shutil.copyfile(outdir + f, jamboreedir + outdir + f)
    


files = []
for (__, __, f) in walk(nonint_dir):
    files += f
for f in files:
    if "nonint" not in f:
        continue
    (patient, sample, gamma, ploidy) = f.split("_")[0:4]
    if onlysomepatients and patient not in somepatients:
        continue
    print("Processing", f)
    nonints = getNonIntegerCalls(f)
    bal_calls = lsl.readBalancedCalls(balance_dir, patient, sample)
    modifyDifficultCalls(nonints, bal_calls)
    writeNewNonInts(f, nonints)
    copyToJamboreeIfWGS(f)

