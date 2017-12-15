# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
#from os.path import isfile
#import matplotlib.pyplot as plt
import lucianSNPLibrary as lsl
import numpy


flist = []
SNPfiles = []
#SNPfiles.append(["1034", "20008", "1034_20008_BAF.txt"])
for (_, _, f) in walk("SNP_all_BAFs/"):
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
    SNPfiles += [[patient, sample, f]]

overlap_files = []
flist = []
for (_, _, f) in walk("CN_calc_log2rs/"):
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
    for snpfile in SNPfiles:
        if (patient == snpfile[0] and sample == snpfile[1]):
            overlap_files += [[f, snpfile[2]]]
            break
        
# read the probeset file, which correlates name to position.
labels, rev_labels = lsl.getSNPLabels()
        
for filepair in overlap_files:
    segfile = filepair[0]
    baffile = filepair[1]
    alldata = {}
    split= segfile.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    bloodfile = patient + "_blood_BAF.txt"
    segments = []
    infile = open("CN_calc_log2rs/" + segfile, "r")
    for line in infile:
        (chr, start, end, x_log2r, call, n_log2r, avg_log2r) = line.rstrip().split()
        if (chr == "chr"):
            continue
        chr = int(chr)
        start = int(start)
        if (end == "inf"):
            end = numpy.inf
        else:
            end = int(end)
        segments.append([chr, start, end])
        
    #Process the BE SNP data
    infile = open("SNP_all_BAFs/" + baffile, "r")
    indata = infile.readlines()
    infile.close()
    labelrow = indata[0].rstrip().split("\t")
    BErow = indata[1].rstrip().split("\t")

    #Process the blood SNP data
    infile = open("SNP_all_BAFs/" + bloodfile, "r")
    indata = infile.readlines()
    infile.close()
    bloodrow = indata[1].rstrip().split("\t")

    SNPs = {}
    for col in range(len(labelrow)):
        if labelrow[col] in labels:
            if (BErow[col] == "?"):
                continue
            if (bloodrow[col] == "?"):
                continue
            SNPs[labels[labelrow[col]]] = (float(BErow[col]), float(bloodrow[col]))

    blood_minus_BE = {}
    #Now go through each SNP, calcuate abs(blood-BE) *if* 0.3<blood<0.7, and save it
    for segment in segments:
        blood_minus_BE[(segment[0], segment[1], segment[2])] = []
        
    for SNP in SNPs.items():
        if (SNP[1][1] >= 0.3 and SNP[1][1] <= 0.7):
            for segment in segments:
                if (SNP[0][0] == segment[0] and SNP[0][1] > segment[1] and SNP[0][1] <= segment[2]):
                    blood_minus_BE[(segment[0], segment[1], segment[2])].append(abs(SNP[1][0] - SNP[1][1]))

    outfile = open("SNP_segment_BAFs/" + patient + "_" + sample + "_segment_BAFdiff.txt", "w")
    outfile.write("chr\tstart\tend\tavg\tmax\tmin\tnumber\n")
    for diff in blood_minus_BE.items():
        outfile.write(str(diff[0][0]) + "\t" + str(diff[0][1]) + "\t" + str(diff[0][2]) + "\t")
        if (len(diff[1]) == 0):
            outfile.write("--\t--\t--\t0\n")
        else:
            outfile.write(str(numpy.average(diff[1])) + "\t" + str(numpy.max(diff[1])) + "\t" + str(numpy.min(diff[1])) + "\t" + str(len(diff[1])) + "\n")
#        for el in diff[1]:
#            outfile.write("\t" + str(el))
#        outfile.write("\n")
