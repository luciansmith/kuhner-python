    # -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 13:57:27 2016

@author: lpsmith
"""

import lucianSNPLibrary as lsl
import string
from os import walk
from os import path
from os import mkdir

# read the file that correlates patient data with which omni file the 'canonical' version of that data can be found.
use_canonical = False
use_averaged_SNPs = False

baf_dir = "BAF_raw_data/"
tag = "_1M"
if use_averaged_SNPs:
    tag += "_averaged"
else:
    tag += "_only"

canonical_filename = "CN_raw_data/20170724_sample_omni.txt"
baf_outdir = "BAF_first_filtered_data" + tag + "/"
baf_labelfile = "REI_12051_B01_SOM_WGS_443samples_12Dec2016_Partek_BAlleleFrequency.1.txt.fmt"

onepatientonly = False
onepatient = "791"

if not(path.isdir(baf_outdir)):
    mkdir(baf_outdir)

# read the probeset file, which correlates name to position.  We dont want zeroes ('False').
if (use_averaged_SNPs):
    labels, rev_labels = lsl.getSNPLabelsAveraged(False)
else:
    labels, rev_labels = lsl.getSNPLabelsAll(False)

if (use_canonical):
    infile = open(canonical_filename,"r")

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
        whichdata[id + "_BLD2"] = experiment
    infile.close()

#filenames = ["omniMix3_BAF.txt"]
filenames = []
for (_, _, f) in walk(baf_dir):
    filenames += f
    break

SNPnames = []
if baf_labelfile in filenames:
    labelfile = open(baf_dir + baf_labelfile, "r")
    start = False
    for line in labelfile:
        if line.find("Sample ID") != -1:
            start = True
        if (line.find("colgroup") != -1):
            continue
        if (start):
            SNPnames.append(line.rstrip())
#print len(SNPnames)
#print filenames

for file in filenames:
    if (file.find("REI") != 0 and file.find("omni") != 0):
        continue
    if file.find(".fmt") != -1:
        continue
    print file
    filebits = file.rstrip().split("_")
    whichomni = filebits[0]
    SNPfile = open(baf_dir + file, "r")
    if baf_labelfile not in filenames:
        SNPnames = SNPfile.readline().rstrip().split("\t")
    for line in SNPfile:
        SNPline = line.rstrip().split("\t")
        id = SNPline[0]
        if (id == "991_23163N_89V_BLD_"):
            id = "991_23163N_89V_BLD"
        idbits = id.rstrip().split("_")
        if (len(idbits) != 4):
            idbits = id.split("-")
        if (len(idbits) != 4):
            print "Not saving data for", SNPline[0], ": label not in canonical form."
            continue
        id = idbits[0]
        if (onepatientonly and id != onepatient):
            continue
        num = idbits[1]
        isblood = idbits[3]
        if (isblood == "BLD"):
            num = "BLD"
        elif (isblood == "gastric"):
            num = "gastric"
        elif (isblood == "BLD2"):
            num = "BLD2"
        elif (isblood == "GASTRIC"):
            num = "gastric"
        elif (isblood == "Gastric"):
            num = "gastric"
        elif (isblood == "GST"):
            num = "gastric"
        if id == "844" and num == "21570":
            num = "21640" #Typo in data entry: 21640 is the correct value; there is no 21570.
#        if id != "991":
#            continue
        idnum = id + "_" + num
#        if (num != "BLD2"):
#            continue
        if (use_canonical):
            if (idnum not in whichdata):
                print "Not saving data for", SNPline[0], ": not in canonical list."
                continue
            refomni = whichdata[idnum]
            if (refomni != whichomni):
                print "The file " + whichomni + " contains data for " + id + "_" + num + ", but that data should be found in " + refomni + "instead."
                continue
        if (len(SNPnames) != len(SNPline)):
            print str(len(SNPnames)) + " SNPnames and " + str(len(SNPline)) + " SNPs: unequal numbers in file " + file + ", for id '" + id + "_" + num + "': skipping."
            continue
        outname = baf_outdir + id + "_" + num + "_BAF.txt"
        if path.isfile(outname):
            print "Already have output for", id, num, ": skipping."
        else:
            print "Writing data for patient", id, ", sample", num
            outfile = open(outname, "w")
            outfile.write("SNPid\tchr\tpos\tBAF\n")
            for entry in range(1,len(SNPnames)):
                label = labels.get(SNPnames[entry])
                if (label == None):
                    continue
                if (label[0] == "UNKNOWN"):
                    continue;
                outfile.write(SNPnames[entry] + "\t" + label[0] + "\t" + label[1] + "\t" + SNPline[entry] + "\n")
            outfile.close()
    SNPfile.close()
