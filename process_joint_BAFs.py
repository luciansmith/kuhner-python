# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:27:59 2016

@author: lpsmith
"""

from __future__ import division
from os import mkdir
from os.path import isfile
import os.path as path
import numpy

# read the processed data from the file Xiaohong sent me, which has segmentation and non-2N bafv values

for subset in ["_25M_v2_"]: #, "_only1M"]
    for ploidy in ["tetraploid", "diploid"]:
        tag = subset + ploidy
    
        joint_file = "joint_seg" + tag + ".txt"
        baf_val_dir= "BAF_filtered_data_25M_15/"
        outdir = "BAF_joint_vals" + tag + "/"
        
        if not(path.isdir(outdir)):
            mkdir(outdir)
        
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
        #    if label != ("862", "18992"):
        #        continue
            outfilename = outdir + patient + "_" + sample + "_avgbafvs.txt"
            if isfile(outfilename):
                print "Already have output for", patient, sample
                continue
            print "Processing patient/sample", patient, sample
            bafv_filename = baf_val_dir + patient + "_" + sample + "_BAF.txt"
            if (not(isfile(bafv_filename))):
                print "Can't find patient/sample", patient, sample
                continue
            bafv_file = open(bafv_filename, "r")
            for line in bafv_file:
                if line.find("chr") != -1:
                    continue
                (id, chr, pos, bafv) = line.strip().split()
                if int(chr) > 22:
                    continue
                pos = int(pos)
                try:
                    bafv = float(bafv)
                except:
                    continue
                if (bafv < 0.5):
                    bafv = 1-bafv
                found = False
                for seg in joint_data[label][chr]:
                    if seg[0] <= pos and seg[1] >= pos:
                        seg[7].append(bafv)
                        found = True
                        break
        #        if not(found):
        #            print patient + "\t" + sample + "\t" + chr + "\t" + str(pos) + "\t" + str(bafv) + "\tnot_found"
        #        else:
        #            print patient + "\t" + sample + "\t" + chr + "\t" + str(pos) + "\t" + str(bafv) + "\tfound"
            outfile = open(outfilename, "w")
            outfile.write("chr\tstartpos\tendpos\trawA\trawB\tintA\tintB\tavgBAF\tnBAF_SNPs\n")
            for chr in range(1,23):
                chr = str(chr)
                for seg in joint_data[label][chr]:
                    if len(seg[7]) > 0:
                        seg[6] = numpy.average(seg[7])
                    else:
                        seg[6] = "NA"
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

