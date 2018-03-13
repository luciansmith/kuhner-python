# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:27:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile
from os import mkdir
import os.path as path
import os
import numpy

#for subset in ["", "_only25M"]: #, "_only1M"]
for subset in ["g250"]: #["25M_v2"]: #, "_only1M"]
    for constraint in ["diploid"]: #, "diploid", "tetraploid"]:

        #tag = "_combined_avSNPs" + subset + "_" + constraint
        tag = "_" + subset + "_" + constraint

        joint_file = "joint_seg" + tag + ".txt"
        CN_25M_val_dir= "CN_filtered_data_25M_only/"
        CN_1M_val_dir= "CN_filtered_data_1M_only/"
        outdir = "CN_joint_log2rs" + tag + "/"

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
            if (rawA == "NA"):
                joint_data[label][chr].append([int(start), int(end), "nan", "nan", "nan", "nan", 0, []])
            else:
                joint_data[label][chr].append([int(start), int(end), float(rawA), float(rawB), int(intA), int(intB), 0, []])


        for label in joint_data:
            (patient, sample) = label
        #    if label != ("521", "18824"):
        #        continue
            if isfile(outdir + patient + "_" + sample + "_avglog2rs.txt"):
                print("Already have output for", patient, sample)
                continue
            print("Processing patient/sample", patient, sample)
            log2r_filename = CN_25M_val_dir + patient + "_" + sample + "_copynumber_all.txt"
            if (not(isfile(log2r_filename))):
                log2r_filename = CN_1M_val_dir + patient + "_" + sample + "_copynumber_all.txt"
                if (not(isfile(log2r_filename))):
                    print("Can't find patient/sample", patient, sample)
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
            filename = outdir + patient + "_" + sample + "_avglog2rs.txt"
            if isfile(filename):
                print("The file", filename, "already exists: assuming it's fine and moving on.")
            else:
                outfile = open(filename, "w")
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

