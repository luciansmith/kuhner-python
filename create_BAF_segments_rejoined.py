#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:02:52 2016

@author: lpsmith
"""

#Create parallel BAF segmentation files to match the CN segmentation files

from __future__ import division
from os import walk
import numpy

import lucianSNPLibrary as lsl

#Use if you want to analyze the 'rejoined' data:
CN_directory = "CN_calc_log2rs_rejoined/"
BAF_directory = "SNP_all_BAFs/"
outDirectory = "BAF_persegment_rejoined/"

#Use if you want to analyze the original Xiaohong-segmented data:
#CN_directory = "CN_calc_log2rs/"
#BAF_directory = "SNP_all_BAFs/"
#outDirectory = "BAF_persegment/"

#Use if you want to analyze the full Xiaohong-segmented data:
#CN_directory = "CN_calc_log2rs/"
#BAF_directory = "SNP_all_BAFs/"
#outDirectory = "BAF_persegment/"

# read the probeset file, which correlates name to position.
labels, rev_labels = lsl.getSNPLabels()

def getBloodData(overlap):
    print "Loading blood data for all samples."
    blooddata = {}
    for filepairs in overlap:
        patient = filepairs[1]
        bloodname = patient + "_blood_BAF.txt"
        if not(bloodname in blooddata):
            hetlist = set()
            bloodfile = open(BAF_directory + bloodname, "r")
            SNPnames = bloodfile.readline().rstrip().split("\t")
            SNPvals = bloodfile.readline().rstrip().split("\t")
            
            for s in range(0, len(SNPnames)):
                try:
                    val = float(SNPvals[s])
                except:
                    continue
                if val >= 0.4 and val <= 0.6:
                    hetlist.add(SNPnames[s])
            blooddata[bloodname] = hetlist
    print "Done loading blood data."
    return blooddata

def getSNPsFor(patient, sample, blooddata):
    print "Getting SNPs for patient", patient, "sample", sample
    ret = {}
    BAFname = patient + "_" + sample + "_BAF.txt"
    bloodname = patient + "_blood_BAF.txt"
    hetlist = blooddata[bloodname]
    bfile = open(BAF_directory + BAFname, "r")
    SNPnames = bfile.readline().rstrip().split("\t")
    SNPvals = bfile.readline().rstrip().split("\t")
    for s in range(0, len(SNPnames)):
        if SNPnames[s] in hetlist and SNPnames[s] in labels:
            label = labels[SNPnames[s]]
            chr = int(label[0])
            pos = int(label[1])
            try:
                baf = float(SNPvals[s])
            except:
                continue
            if not(chr in ret):
                ret[chr] = []
            ret[chr].append([pos, baf])
    print "Obtained SNPs."
    return ret

flist = []
SNPfiles = []
#SNPfiles.append(["1034", "20008"])
for (_, _, f) in walk(BAF_directory):
    flist += f

for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (sample == "blood"):
        continue
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    SNPfiles += [[patient, sample]]

doubled = [[141, 21060], [141, 21062], [141, 21064], [163, 19208], [163, 19214], [194, 19868], [194, 19880], [450, 18974], [450, 18982], [512, 18744], [512, 18746], [512, 18748], [512, 18750], [512, 18762], [660, 19260], [660, 19262], [660, 19264], [660, 19272], [664, 19954], [772, 18944], [772, 18946], [848, 18794], [884, 20354], [884, 20358], [954, 20014], [954, 20016], [954, 20018], [991, 20600], [997, 20656], [997, 20666], [997, 20668], [997, 20672], [997, 20674], [1006, 21104], [1044, 20856], [1044, 20864], [997, 20658], [997, 20660], [660, 19266], [660, 19270], [740, 20000], [997, 20662], [997, 20664]]


overlap_files = []
flist = []
for (_, _, f) in walk(CN_directory):
    flist += f

for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    if [int(patient), int(sample)] in doubled:
        continue
    if [patient, sample] in SNPfiles:
        overlap_files += [[f, patient, sample]]


#Load the blood data
blooddata = getBloodData(overlap_files)
#print "blooddata = ", blooddata

for filepairs in overlap_files:
    patient = int(filepairs[1])
    sample = int(filepairs[2])
    if (sample != 20008):
        continue
    CNfile = open(CN_directory + filepairs[0], "r")
    bafdata = getSNPsFor(filepairs[1], filepairs[2], blooddata)
    #print "bafdata: ", bafdata
    total_n = 0
    segments = {}
    for line in CNfile:
        if (line.find("chr") != -1):
            continue
        (chr, start, end, x_log2r, call, n_log2r, avg_log2r, stdev) = line.rstrip().split()
        chr = int(chr)
        if (end=="inf"):
            end = float('inf')
        else:
            end = int(end)
        if not(chr in segments):
            segments[chr] = []
        segments[chr].append([int(start), end, avg_log2r, list()])
    for bchr in range(1, 23):
        for (bpos, val) in bafdata[bchr]:
            for segment in segments[bchr]:
                if bpos <= segment[0]:
                    continue
                if bpos > segment[1]:
                    continue
                if (val < 0.5):
                    val = 1-val
                segment[3].append(val)
                break;
    print "Writing output for patient", patient, "sample", sample
    outfile = open(outDirectory + filepairs[1] + "_" + filepairs[2] + ".txt", "w")
    outfile.write("chr\tstart\tend\tavg_log2r\tavg_BAF\n")
    for c in range(1,23):
        for segment in segments[c]:
            outfile.write(str(c) + "\t")
            outfile.write(str(segment[0]) + "\t")
            outfile.write(str(segment[1]) + "\t")
            outfile.write(segment[2] + "\t")
            if len(segment[3]) > 0:
                outfile.write(str(numpy.average(segment[3])))
                if (c == 3 and segment[0] == 60384729):
                    print "BAF values for chr3,", segment[0], ",", segment[1]
                    print segment[3]
            else:
                outfile.write("---")
            outfile.write("\n")
    outfile.close()



