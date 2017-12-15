# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 13:57:27 2016

@author: lpsmith
"""

import lucianSNPLibrary as lsl
from os import walk
import string

# read the probeset file, which correlates name to position.
labels, rev_labes = lsl.getSNPLabels()

# read the file that correlates patient data with which omni file the 'canonical' version of that data can be found.
infilename = "BAF_raw_data/20160621_sample_omni.txt"
infile = open(infilename,"r")

whichdata = {}
for line in infile:
    if (line[0]=="PatientID"):
        continue
    line = line.rstrip().split("\t")
    (id, num, experiment) = line
    whichdata[id + "_" + num] = experiment
    if id + "_BLD" in experiment:
        if whichdata[id + "_BLD"] != experiment:
            print id + "_BLD potentially in two different files:", experiment, "and", whichdata[id + "_BLD"]
    whichdata[id + "_BLD"] = experiment
    whichdata[id + "_gastric"] = experiment
infile.close()

#filenames = ["omni15_paired_copynumber_baselineAdjusted_37_RB2012_2015.txt"]
filenames = []
for (_, _, f) in walk("BAF_raw_data/"):
    filenames += f
    break

for file in filenames:
    if (string.find(file, "omni") != 0):
        continue
    filebits = file.rstrip().split("_")
    whichomni = filebits[0]
    SNPfile = open("BAF_raw_data/" + file,"r")
    SNPnames = SNPfile.readline().rstrip().split("\t")
    for line in SNPfile:
        SNPline = line.rstrip().split("\t")
        id = SNPline[0]
        idbits = id.rstrip().split("_")
        if (len(idbits) != 4):
            print SNPline[0] + "\t" + file + "\tNot canonical"
            continue
        id = idbits[0]
        num = idbits[1]
        isblood = idbits[3]
        if (isblood == "BLD"):
            num = "BLD"
        if (isblood == "BLD2"):
            num = "BLD"
        elif (isblood == "gastric"):
            num = "gastric"
        elif (isblood == "GASTRIC"):
            num = "gastric"
        elif (isblood == "Gastric"):
            num = "gastric"
        idnum = id + "_" + num
        if (idnum not in whichdata):
            print SNPline[0] + "\t" + file + "\tNot listed"
            continue
        refomni = whichdata[idnum]
        if (refomni != whichomni):
            print SNPline[0] + "\t" + file + "\tRe-ran"
            continue
        if (len(SNPnames) != len(SNPline)):
            print SNPline[0] + "\t" + file + "\tProblem with file"
            continue
        print SNPline[0] + "\t" + file + "\tListed"
        whichdata[idnum] = "found"

for idnum in whichdata:
    if whichdata[idnum] != "found":
        print idnum + "\t" + whichdata[idnum] + "\tNot found"
