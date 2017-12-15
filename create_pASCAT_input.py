#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:06:34 2017

@author: lpsmith
"""

import lucianSNPLibrary as lsl
from os import walk
from os import path
from os import mkdir

baf_dir = "BAF_first_filtered_data_25M/"
cn_dir = "CN_filtered_data_25M/"
outdir = "pASCAT_input_25M/"

if not(path.isdir(outdir)):
    mkdir(outdir)

onepatientonly = False
onepatient = "450"


flist = []
CNfiles = {}
for (_, _, f) in walk(cn_dir):
    flist += f
   
for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if onepatientonly and patient!=onepatient:
        continue
    if not patient in CNfiles:
        CNfiles[patient] = {}
    CNfiles[patient][sample] = f


flist = []
baffiles = {}
for (_, _, f) in walk(baf_dir):
    flist += f
   
for f in flist:
    if (f.find(".txt") == -1):
        continue
    split = f.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if onepatientonly and patient!=onepatient:
        continue
    if not patient in baffiles:
        baffiles[patient] = {}
    baffiles[patient][sample] = f

removedpatients = []
removedsamples = []
for patient in CNfiles:
    if not patient in baffiles:
        removedpatients.append(patient)
        print "Removing", patient, "from CNfiles: no such patient in BAF files."
        continue
    for sample in CNfiles[patient]:
        if not sample in baffiles[patient]:
            print "Removing", patient, ",", sample, "from CNfiles: no such sample in BAF files."
            removedsamples.append((patient, sample))
    

for patient in baffiles:
    if not patient in CNfiles:
        removedpatients.append(patient)
        print "Removing", patient, "from baffiles: no such patient in BAF files."
        continue
    for sample in baffiles[patient]:
        if (sample == "gastric" or sample=="BLD" or sample.find("N") != -1):
            continue
        if not sample in CNfiles[patient]:
            print "Removing", patient, ",", sample, "from baffiles: no such sample in CN files."
            removedsamples.append((patient, sample))
            
for patient in removedpatients:
    if patient in CNfiles:
        del CNfiles[patient]
    if patient in baffiles:
        del baffiles[patient]

for (patient, sample) in removedsamples:
    if sample in CNfiles[patient]:
        del CNfiles[patient][sample]
    if sample in baffiles[patient]:
        del baffiles[patient][sample]

for patient in CNfiles:
    if onepatientonly and patient != onepatient:
        continue
    CN_data = {}
    samples = []
    print "Processing CN output for", patient
    for sample in CNfiles[patient]:
        cnfilename = CNfiles[patient][sample]
        if (cnfilename.find("N") != -1):
            print "Skipping", cnfilename, ": probably a gastric sample."
            continue
        samples.append(sample)
        cnfile = open(cn_dir + cnfilename, "r")
        print "Reading", cnfilename
        for line in cnfile:
            if "SNPid" in line:
                continue
            (id, chr, pos, l2r) = line.rstrip().split("\t")
            try:
                float(l2r)
            except:
                if l2r != "?":
                    print "Non-float value", l2r, "for", id
                l2r = "?"
            if chr == "23" or chr=="24":
                continue
            id = lsl.fixNameForR(id)
            if not id in CN_data:
                CN_data[id] = [chr, pos, {}]
            CN_data[id][2][sample] = l2r
    cnout = open(outdir + patient + "_logR.txt", "w")
    cnout.write("\t\"Chr\"\t\"Position\"")
    for sample in samples:
        cnout.write("\t\"" + patient + "_" + sample + "\"")
    cnout.write("\n")
    sorted_CNs = []
    for id in CN_data:
        (chr, pos, sset) = CN_data[id]
        sorted_CNs.append((int(chr), int(pos), id, sset))
    sorted_CNs.sort()
    for (chr, pos, id, sset) in sorted_CNs:
        if id=="rs2098322":
            continue
#    for id in CN_data:
#        (chr, pos, sset) = CN_data[id]
        cnout.write("\"" + id + "\"")
        cnout.write("\t\"" + str(chr) + "\"\t" + str(pos))
        for sample in samples:
            l2r = sset[sample]
            if (l2r == "?"):
                cnout.write("\tNA")
            else:
                cnout.write("\t" + str(float(l2r)))
        cnout.write("\n")
    cnout.close()

for patient in baffiles:
    if onepatientonly and patient != onepatient:
        continue
    baf_data = {}
    samples = []
    bloodorgastric = ""
    print "Processing baf output for", patient
    for sample in baffiles[patient]:
        baffilename = baffiles[patient][sample]
        if (sample == "gastric" or sample=="BLD" or sample.find("N") != -1):
            bloodorgastric = sample
        else:
            samples.append(sample)
        baffile = open(baf_dir + baffilename, "r")
        print "Reading", baffilename
        for line in baffile:
            if "SNPid" in line:
                continue
            (id, chr, pos, l2r) = line.rstrip().split("\t")
            if chr == "23" or chr=="24":
                continue
            id = lsl.fixNameForR(id)
            if not id in baf_data:
                baf_data[id] = [chr, pos, {}]
            baf_data[id][2][sample] = l2r
    bafout = open(outdir + patient + "_BAF.txt", "w")
    bafNout = open(outdir + patient + "_Normal_BAF.txt", "w")
    lineout = "\t\"Chr\"\t\"Position\""
    for sample in samples:
        lineout += "\t\"" + patient + "_" + sample + "\""
    lineout += "\n"
    bafout.write(lineout)
    bafNout.write(lineout)
    sorted_bafs = []
    for id in baf_data:
        (chr, pos, sset) = baf_data[id]
        sorted_bafs.append((int(chr), int(pos), id, sset))
    sorted_bafs.sort()
    for (chr, pos, id, sset) in sorted_bafs:
        if id=="rs2098322":
            continue
        bafout.write("\"" + id + "\"\t\"" + str(chr) + "\"\t" + str(pos))
        bafNout.write("\"" + id + "\"\t\"" + str(chr) + "\"\t" + str(pos))
        for sample in samples:
            l2r = sset[sample]
            if (l2r == "?"):
                bafout.write("\tNA")
            else:
                bafout.write("\t" + str(float(l2r)))
            l2r = sset[bloodorgastric]
            if (l2r == "?"):
                bafNout.write("\tNA")
            else:
                bafNout.write("\t" + str(float(l2r)))
        bafout.write("\n")
        bafNout.write("\n")
    bafout.close()
    bafNout.close()

