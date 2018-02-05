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

bafdirs = {}
cndirs = {}


use_averaged_SNPs = False

onepatientonly = False
onepatient = "572"

if use_averaged_SNPs:
    labels, rev_labels = lsl.getSNPLabelsAveraged(False)
    avtag = "avSNPs"
    dirtag = "_averaged"
else:
    labels, rev_labels = lsl.getSNPLabelsAll(False)
    avtag = "all"
    dirtag = "_only"

bafdirs["25"] = "BAF_first_filtered_data_25M" + dirtag + "/"
cndirs["25"] = "CN_filtered_data_25M" + dirtag + "/"
bafdirs["1"] = "BAF_first_filtered_data_1M" + dirtag + "/"
cndirs["1"] = "CN_filtered_data_1M" + dirtag + "/"


def readBAFFile(baffilename, sample, avgstraight):
    baffile = open(baffilename, "r")
    print "Reading", baffilename
    for line in baffile:
        if "SNPid" in line:
            continue
        (id, chr, pos, l2r) = line.rstrip().split("\t")
        if id.find("cnvi") != -1:
            continue
        if use_averaged_SNPs:
            if id not in labels:
                continue
            #The SNP might have changed locations
            (chr, pos) = labels[id]
        if chr == "23" or chr=="24":
            continue
        if not chr in baf_data:
            baf_data[chr] = {}
        if not pos in baf_data[chr]:
            baf_data[chr][pos] = {}
        if sample in baf_data[chr][pos]:
            #multiple data points for same sample.  Assuming two, and averaging:
            if l2r == "?":
                #Don't save the data; leave it as-is.
                continue
            if baf_data[chr][pos][sample] == "?":
                baf_data[chr][pos][sample] = l2r
            else:
                l2r = float(l2r)
                avgwith = float(baf_data[chr][pos][sample])
                #print "Averaging baf data for", patient, sample, chr, pos, l2r, avgwith
                if not(avgstraight) and l2r<0.5 and avgwith>0.5:
                    l2r = 1-l2r
                elif not(avgstraight) and l2r > 0.5 and avgwith<0.5:
                    avgwith = 1-avgwith
                l2r = str((l2r+avgwith)/2)
        baf_data[chr][pos][sample] = l2r

for (only25, only1) in [(False, False), (False, True), (True, False)]:
    onlytag = ""
    if (only25):
        onlytag = "_only25M"
    if (only1):
        onlytag = "_only1M"

    outdir = "pASCAT_input_combined_" + avtag + onlytag + "/"
    print "Creating output in", outdir
    if not(path.isdir(outdir)):
        mkdir(outdir)

    dirs = ["1", "25"]
    if (only25):
        dirs = ["25"]
    elif (only1):
        dirs = ["1"]


    CNfiles = {}
    baffiles = {}

    patients = {}
    patients["25"] = set()
    patients["1"] = set()

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
            if onepatientonly and patient!=onepatient:
                continue
            if not patient in CNfiles:
                CNfiles[patient] = {}
            if sample not in CNfiles[patient]:
                CNfiles[patient][sample] = CNdir + f
            else:
                CNfiles[patient][sample] = [CNfiles[patient][sample], CNdir + f]

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
            if onepatientonly and patient!=onepatient:
                continue
            if not patient in baffiles:
                baffiles[patient] = {}
            if sample == "BLD" or sample=="BLD2" or sample=="gastric" or sample=="Blood" or sample.find("N") != -1:
                sample = "Blood"
                if sample not in baffiles[patient]:
                    baffiles[patient][sample] = []
                baffiles[patient][sample].append(bafdir + f)
            elif sample not in baffiles[patient]:
                baffiles[patient][sample] = bafdir + f
            else:
                print "Multiple files for", patient, sample
                baffiles[patient][sample] = [baffiles[patient][sample], bafdir + f]

    removedpatients = []
    removedsamples = []

    if not only25 and not only1:
        # We're combining the files; only keep patients that have data for both
        for patient in patients["25"]:
            if patient not in patients["1"]:
                print "Removing patient", patient, "from combined analysis: only have samples with 2.5M data"
                removedpatients.append(patient)
        for patient in patients["1"]:
            if patient not in patients["25"]:
                print "Removing patient", patient, "from combined analysis: only have samples with 1M data"
                removedpatients.append(patient)


    for patient in CNfiles:
        if not patient in baffiles:
            removedpatients.append(patient)
            print "Removing", patient, "from CNfiles: no such patient in BAF files."
            continue
        if len(CNfiles[patient].keys()) < 3:
            removedpatients.append(patient)
            print "Removing", patient, "from CNfiles: too few samples (<3)."
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
        if len(baffiles[patient].keys()) < 4:
            removedpatients.append(patient)
            print "Removing", patient, "from baffiles: too few samples (<3)."
            continue
        for sample in baffiles[patient]:
            if (sample == "gastric" or sample=="BLD" or sample=="BLD2" or sample=="Blood" or sample.find("N") != -1):
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
        if patient in removedpatients:
            continue
        if sample in CNfiles[patient]:
            del CNfiles[patient][sample]
        if sample in baffiles[patient]:
            del baffiles[patient][sample]

    for patient in CNfiles:
        if onepatientonly and patient != onepatient:
            continue
        CN_data = {}
        samples = []
        if path.isfile(outdir + patient + "_logR.txt"):
            print "Skipping", patient, "logR: file already exists"
            continue
        print "Processing CN output for", patient
        for sample in CNfiles[patient]:
            cnfilename = CNfiles[patient][sample]
            if (cnfilename.find("N_c") != -1):
                print "Skipping", cnfilename, ": probably a gastric sample."
                continue
            samples.append(sample)
            cnfile = open(cnfilename, "r")
            print "Reading", cnfilename
            for line in cnfile:
                if "SNPid" in line:
                    continue
                (id, chr, pos, l2r) = line.rstrip().split("\t")
                if chr == "23" or chr=="24":
                    continue
                if use_averaged_SNPs:
                    #The SNP might have changed locations
                    if id not in labels:
                        continue
                    (chr, pos) = labels[id]
                chr = int(chr)
                pos = int(pos)
                try:
                    float(l2r)
                except:
                    if l2r != "?":
                        print "Non-float value", l2r, "for", id
                    l2r = "?"
                if chr == "23" or chr=="24":
                    continue
                if not chr in CN_data:
                    CN_data[chr] = {}
                if not pos in CN_data[chr]:
                    CN_data[chr][pos] = {}
                if sample in CN_data[chr][pos]:
                    #multiple data points for same sample.  Assuming two, and averaging:
                    if l2r == "?":
                        #Don't save the data; leave it as-is.
                        continue
                    if CN_data[chr][pos][sample] == "?":
                        CN_data[chr][pos][sample] = l2r
                    else:
    #                    print "Averaging CN data for", patient, sample, chr, pos
                        l2r = str((float(l2r) + float(CN_data[chr][pos][sample]))/2)
                CN_data[chr][pos][sample] = l2r
        cnout = open(outdir + patient + "_logR.txt", "w")
        cnout.write("\t\"Chr\"\t\"Position\"")
        for sample in samples:
            cnout.write("\t\"" + patient + "_" + sample + "\"")
        cnout.write("\n")
        chr_keys = CN_data.keys()
        chr_keys.sort()
        for chr in chr_keys:
            pos_keys = CN_data[chr].keys()
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
        if onepatientonly and patient != onepatient:
            continue
        baf_data = {}
        samples = []
        bloodorgastric = "Blood"
        if path.isfile(outdir + patient + "_BAF.txt"):
            print "Skipping", patient, "BAF: file already exists"
            continue
        print "Processing baf output for", patient
        for sample in baffiles[patient]:
            if sample == "Blood":
                for filename in baffiles[patient][sample]:
                    print "One Blood sample found", filename
                    readBAFFile(filename, sample, True)
            else:
                baffilename = baffiles[patient][sample]
                samples.append(sample)
                readBAFFile(baffilename, sample, False)


        bafout = open(outdir + patient + "_BAF.txt", "w")
        bafNout = open(outdir + patient + "_Normal_BAF.txt", "w")
        lineout = "\t\"Chr\"\t\"Position\""
        for sample in samples:
            lineout += "\t\"" + patient + "_" + sample + "\""
        lineout += "\n"
        bafout.write(lineout)
        bafNout.write(lineout)
        chr_keys = baf_data.keys()
        chr_keys.sort()
        for chr in chr_keys:
            pos_keys = baf_data[chr].keys()
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

