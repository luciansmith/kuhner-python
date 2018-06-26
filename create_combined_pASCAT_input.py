#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:06:34 2017

@author: lpsmith
"""

from __future__ import division
import lucianSNPLibrary as lsl
from os import walk
from os import path
from os import mkdir
import numpy


somepatientsonly = False
somepatients = ["997"]
#somepatients = ["551", "59", "954", "222", "88", "639", "611", "391", "422", "575", "619", "672", "728", "915", "1005", "43", "686"]

labels, rev_labels = lsl.getSNPLabelsAll(False)
avtag = "all"
dirtag = "_only"

bafdirs = {}
cndirs = {}
bafdirs["25"] = "BAF_first_filtered_data_25M" + dirtag + "/"
cndirs["25"] = "CN_filtered_data_25M" + dirtag + "/"
bafdirs["1"] = "BAF_first_filtered_data_1M" + dirtag + "/"
cndirs["1"] = "CN_filtered_data_1M" + dirtag + "/"
bafdirs["Pilot"] = "BAF_first_filtered_data_Pilot/"
cndirs["Pilot"] = "CN_filtered_data_Pilot/"

#bad_dupes = open("bad_duplicates.txt", "w")
#bad_dupes.write("Patient\tSample\tChr\tpos\tval1\tval2\n")

#blood_diffs = []

def readBAFFile(baffilename, sample, isblood, baf_data):
    baffile = open(baffilename, "r")
    print("Reading", baffilename)
    for line in baffile:
        if "SNPid" in line:
            continue
        if "cnvi" in line:
            continue
        (id, chr, pos, bafval) = line.rstrip().split("\t")
        if chr == "23" or chr=="24":
            continue
        if bafval == "?":
            continue
        bafval = float(bafval)
        if not chr in baf_data:
            baf_data[chr] = {}
        pos = int(pos)
        if not pos in baf_data[chr]:
            baf_data[chr][pos] = {}
        if sample not in baf_data[chr][pos]:
            baf_data[chr][pos][sample] = []
        baf_data[chr][pos][sample].append(bafval)

def readCNFile(cnfilename, sample, CN_data):
    cnfile = open(cnfilename, "r")
    print("Reading", cnfilename)
    for line in cnfile:
        if "SNPid" in line:
            continue
        if "cnvi" in line:
            continue
        (id, chr, pos, l2r) = line.rstrip().split("\t")
        if chr == "23" or chr=="24":
            continue
        if l2r == "?":
            continue
        chr = int(chr)
        pos = int(pos)
        try:
            l2r = float(l2r)
            if cnfilename == bafdirs["Pilot"] + "391_19578_copynumber_all.txt":
                l2r = float(l2r) - 0.32260464413499457
        except:
            print("Non-float value", l2r, "for", id)
            assert(False)
            continue
        if not chr in CN_data:
            CN_data[chr] = {}
        if not pos in CN_data[chr]:
            CN_data[chr][pos] = {}
        if sample not in CN_data[chr][pos]:
            CN_data[chr][pos][sample] = []
        CN_data[chr][pos][sample].append(l2r)

def averageData(data, isbaf):
    for chr in data:
        for pos in data[chr]:
            for sample in data[chr][pos]:
                vec = data[chr][pos][sample]
                too_different = False
                if isbaf and len(vec)>1:
                    for n in range(1,len(vec)):
                        if (vec[n] > 0.5 and vec[0] < 0.5) or (vec[n] < 0.5 and vec[0] > 0.5):
                            vec[n] = 1-vec[n]
                    for n in range(0,len(vec)):
                        for m in range(n, len(vec)):
                            if abs(vec[n] - vec[m]) > 0.35:
                                too_different = True
                if too_different:
                    data[chr][pos][sample] = "?"
                else:
                    if len(vec)>2:
                        print(vec)
                    data[chr][pos][sample] = numpy.average(vec)

for (only25, only1) in [(False, False)]:#, (False, True), (True, False)]:
    onlytag = ""
    if (only25):
        onlytag = "_only25M"
    if (only1):
        onlytag = "_only1M"

    outdir = "pASCAT_input_combined_" + avtag + onlytag + "/"
    print("Creating output in", outdir)
    if not(path.isdir(outdir)):
        mkdir(outdir)

    dirs = ["1", "25", "Pilot"]
    if (only25):
        dirs = ["25"]
    elif (only1):
        dirs = ["1"]


    CNfiles = {}
    baffiles = {}

    patients = {}
    patients["25"] = set()
    patients["1"] = set()
    patients["Pilot"] = set()

    for onedir in dirs:
        CNdir = cndirs[onedir]
        flist = []
        for (_, _, f) in walk(CNdir):
            flist += f

        for f in flist:
            if (f.find(".txt") == -1):
                continue
            split = f.split("_")
            if (len(split) < 3):
                continue
            patient = split[0]
            sample = split[1]
            patients[onedir].add(patient)
            if somepatientsonly and patient not in somepatients:
                continue
            if not patient in CNfiles:
                CNfiles[patient] = {}
            if sample not in CNfiles[patient]:
                CNfiles[patient][sample] = [CNdir + f]
            else:
                print("Found two infiles for", patient, sample)
                CNfiles[patient][sample].append(CNdir + f)

    for onedir in dirs:
        bafdir = bafdirs[onedir]
        flist = []
        for (_, _, f) in walk(bafdir):
            flist += f

        for f in flist:
            if (f.find(".txt") == -1):
                continue
            split = f.split("_")
            if (len(split) < 3):
                continue
            patient = split[0]
            sample = split[1]
            if somepatientsonly and patient not in somepatients:
                continue
            if not patient in baffiles:
                baffiles[patient] = {}
            if sample == "BLD" or sample=="BLD2" or sample=="gastric" or sample=="Blood":
                sample = "Blood"
                if sample not in baffiles[patient]:
                    baffiles[patient][sample] = []
                baffiles[patient][sample].append(bafdir + f)
            elif sample not in baffiles[patient]:
                baffiles[patient][sample] = [bafdir + f]
            else:
                print("Multiple files for", patient, sample)
                baffiles[patient][sample].append(bafdir + f)

    removedpatients = []
    removedsamples = []

#    if not only25 and not only1:
#        # We're combining the files; only keep patients that have data for both
#        for patient in patients["25"]:
#            if patient not in patients["1"]:
#                print("Removing patient", patient, "from combined analysis: only have samples with 2.5M data")
#                removedpatients.append(patient)
#        for patient in patients["1"]:
#            if patient not in patients["25"]:
#                print("Removing patient", patient, "from combined analysis: only have samples with 1M data")
#                removedpatients.append(patient)


    for patient in CNfiles:
        if not patient in baffiles:
            removedpatients.append(patient)
            print("Removing", patient, "from CNfiles: no such patient in BAF files.")
            continue
        if len(CNfiles[patient].keys()) < 3:
            removedpatients.append(patient)
            print("Removing", patient, "from CNfiles: too few samples (<3).")
            continue
        for sample in CNfiles[patient]:
            if not sample in baffiles[patient]:
                print("Removing", patient, ",", sample, "from CNfiles: no such sample in BAF files.")
                removedsamples.append((patient, sample))


    for patient in baffiles:
        if not patient in CNfiles:
            removedpatients.append(patient)
            print("Removing", patient, "from baffiles: no such patient in BAF files.")
            continue
        if len(baffiles[patient].keys()) < 4:
            removedpatients.append(patient)
            print("Removing", patient, "from baffiles: too few samples (<3).")
            continue
        for sample in baffiles[patient]:
            if (sample == "gastric" or sample=="BLD" or sample=="BLD2" or sample=="Blood"):
                continue
            if not sample in CNfiles[patient]:
                print("Removing", patient, ",", sample, "from baffiles: no such sample in CN files.")
                removedsamples.append((patient, sample))

    for patient in removedpatients:
        if patient in CNfiles:
            del CNfiles[patient]
        if patient in baffiles:
            del baffiles[patient]

    for (patient, sample) in removedsamples:
        if patient in removedpatients:
            continue
        if sample in CNfiles[patient]:
            del CNfiles[patient][sample]
        if sample in baffiles[patient]:
            del baffiles[patient][sample]

    for patient in CNfiles:
        if somepatientsonly and patient not in somepatients:
            continue
        CN_data = {}
        samples = []
        if path.isfile(outdir + patient + "_logR.txt"):
            print("Skipping", patient, "logR: file already exists")
            continue
        print("Processing CN output for", patient)
        for sample in CNfiles[patient]:
            cnfilenames = CNfiles[patient][sample]
#            if (cnfilename.find("N_c") != -1):
#                print("Skipping", cnfilename, ": probably a gastric sample.")
#                continue
            samples.append(sample)
            for cnfilename in cnfilenames:
                readCNFile(cnfilename, sample, CN_data)
        averageData(CN_data, False)
        cnout = open(outdir + patient + "_logR.txt", "w")
        cnout.write("\t\"Chr\"\t\"Position\"")
        for sample in samples:
            cnout.write("\t\"" + patient + "_" + sample + "\"")
        cnout.write("\n")
        chr_keys = list(CN_data.keys())
        chr_keys.sort()
        for chr in chr_keys:
            pos_keys = list(CN_data[chr].keys())
            pos_keys.sort()
            for pos in pos_keys:
                chrpos = (str(chr), str(pos))
                if chrpos not in rev_labels:
                    #It used to be a 1M label that was given a '0' in 2.5M
                    continue
                id = rev_labels[chrpos]
    #            if id=="rs2098322":
    #                continue
                somedata = False
                linestr = "\"" + id + "\"\t\"" + str(chr) + "\"\t" + str(pos)
                for sample in samples:
                    l2r = "?"
                    if sample in CN_data[chr][pos]:
                        l2r = CN_data[chr][pos][sample]
                    if (l2r == "?"):
                        linestr += "\tNA"
                    else:
                        somedata = True
                        linestr += "\t" + str(float(l2r))
                linestr += "\n"
                if somedata:
                    cnout.write(linestr)
        cnout.close()

    for patient in baffiles:
        if somepatientsonly and patient not in somepatients:
            continue
        baf_data = {}
        samples = []
        bloodorgastric = "Blood"
        if path.isfile(outdir + patient + "_BAF.txt"):
            print("Skipping", patient, "BAF: file already exists")
            continue
        print("Processing baf output for", patient)
        for sample in baffiles[patient]:
            if sample == "Blood":
                for filename in baffiles[patient][sample]:
                    print("One Blood sample found", filename)
                    readBAFFile(filename, sample, True, baf_data)
            else:
                baffilenames = baffiles[patient][sample]
                samples.append(sample)
                for baffilename in baffilenames:
                    readBAFFile(baffilename, sample, False, baf_data)

        averageData(baf_data, True)
        bafout = open(outdir + patient + "_BAF.txt", "w")
        bafNout = open(outdir + patient + "_Normal_BAF.txt", "w")
        lineout = "\t\"Chr\"\t\"Position\""
        for sample in samples:
            lineout += "\t\"" + patient + "_" + sample + "\""
        lineout += "\n"
        bafout.write(lineout)
        bafNout.write(lineout)
        chr_keys = list(baf_data.keys())
        chr_keys.sort()
        for chr in chr_keys:
            pos_keys = list(baf_data[chr].keys())
            pos_keys.sort()
            for pos in pos_keys:
                chrpos = (str(chr), str(pos))
                if chrpos not in rev_labels:
                    #It used to be a 1M label that was given a '0' in 2.5M
                    continue
                id = rev_labels[chrpos]
    #            if id=="rs2098322":
    #                continue
                somedata = False
                linestr = "\"" + id + "\"\t\"" + str(chr) + "\"\t" + str(pos)
                linestrN = "\"" + id + "\"\t\"" + str(chr) + "\"\t" + str(pos)
                for sample in samples:
                    l2r = "?"
                    if sample in baf_data[chr][pos]:
                        l2r = baf_data[chr][pos][sample]
                    if (l2r == "?"):
                        linestr += "\tNA"
                    else:
                        somedata = True
                        linestr += "\t" + str(float(l2r))
                    l2r = "?"
                    if bloodorgastric in baf_data[chr][pos]:
                        l2r = baf_data[chr][pos][bloodorgastric]
                    if (l2r == "?"):
                        linestrN += "\tNA"
                    else:
                        linestrN += "\t" + str(float(l2r))
                if somedata:
                    bafout.write(linestr + "\n")
                    bafNout.write(linestrN + "\n")
        bafout.close()
        bafNout.close()

