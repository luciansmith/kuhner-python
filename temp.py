# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:27:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile
import numpy

# read the processed data from the file Xiaohong sent me, which has segmentation and non-2N log2R values

tag = "_v5_diploid"

joint_file = "joint_seg" + tag + ".txt"
CN_val_dir= "CN_filtered_data/"
outdir = "CN_joint_log2rs" + tag + "/"

joint_data = {}
jfile = open(joint_file, "r")
for line in jfile:
    if line.find("patient") != -1:
        continue
    (patient, sample, chr, start, end, rawA, rawB, intA, intB) = line.strip().split()
    label = (patient, sample)
    if not(label in joint_data):
        joint_data[label] = {}
    if not(chr in joint_data[label]):
        joint_data[label][chr] = []
    joint_data[label][chr].append([int(start), int(end), float(rawA), float(rawB), int(intA), int(intB), 0, []])


for label in joint_data:
    (patient, sample) = label
#    if label != ("521", "18824"):
#        continue
    print "Processing patient/sample", patient, sample
    log2r_filename = CN_val_dir + patient + "_" + sample + "_copynumber_all.txt"
    if (not(isfile(log2r_filename))):
        print "Can't find patient/sample", patient, sample
        continue
    log2r_file = open(log2r_filename, "r")
    for line in log2r_file:
        if line.find("chr") != -1:
            continue
        (id, chr, pos, log2r) = line.strip().split()
        if int(chr) > 22:
            continue
        pos = int(pos)
        try:
            log2r = float(log2r)
        except:
            continue
        
        found = False
        for seg in joint_data[label][chr]:
            if seg[0] <= pos and seg[1] >= pos:
                seg[7].append(log2r)
                found = True
                break
#        if not(found):
#            print patient + "\t" + sample + "\t" + chr + "\t" + str(pos) + "\t" + str(log2r) + "\tnot_found"
#        else:
#            print patient + "\t" + sample + "\t" + chr + "\t" + str(pos) + "\t" + str(log2r) + "\tfound"
    outfile = open(outdir + patient + "_" + sample + "_avglog2rs.txt", "w")
    outfile.write("chr\tstartpos\tendpos\trawA\trawB\tintA\tintB\tavglog2r\tnSNPs\n")
    for chr in range(1,23):
        chr = str(chr)
        for seg in joint_data[label][chr]:
            seg[6] = numpy.average(seg[7])
            outfile.write(chr + "\t")
            outfile.write(str(seg[0]) + "\t")
            outfile.write(str(seg[1]) + "\t")
            outfile.write(str(seg[2]) + "\t")
            outfile.write(str(seg[3]) + "\t")
            outfile.write(str(seg[4]) + "\t")
            outfile.write(str(seg[5]) + "\t")
            outfile.write(str(seg[6]) + "\t")
            outfile.write(str(len(seg[7])) + "\n")
    outfile.close()
        
