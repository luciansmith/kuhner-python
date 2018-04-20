# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 13:57:27 2016

@author: lpsmith
"""

import lucianSNPLibrary as lsl
from os import walk
from os import path
from os import mkdir

# read the file that correlates patient data with which omni file the 'canonical' version of that data can be found.
use_canonical = False
checkonly = False
somepatientsonly = True
somepatients = ["391","611"]
#somepatients = ["422", "575", "619", "672", "728", "915", "1005", "43", "686", "611"]

CN_raw_dir = "CN_raw_data_Pilot/"
tag = "_Pilot"

canonical_filename = "CN_raw_data/20170724_sample_omni.txt"
CN_out_dir = "CN_filtered_data" + tag + "/"

if not(path.isdir(CN_out_dir)):
    mkdir(CN_out_dir)

# read the probeset file, which correlates name to position.  We dont want zeroes ('False').
labels, rev_labels = lsl.getSNPLabelsAll(False)

nan_probes = {}
for probefile in ["nan_probes.txt", "nan_probes2.txt", "nan_probes3.txt"]:
    if (path.isfile(CN_raw_dir + probefile)):
        nanfile = open(CN_raw_dir + probefile)
        for line in nanfile:
            if line.find("Batch") != -1:
                continue
            linevec = line.rstrip().split("\t")
            sample = linevec[1].split("-")[1]
            nan_probes[sample] = linevec[2:]

infile = open(canonical_filename,"r")

whichdata = {}
if (use_canonical):
    for line in infile:
        if (line[0]=="PatientID"):
            continue
        line = line.rstrip().split("\t")
        (id, sample, experiment) = line
        whichdata[id + "_" + sample] = experiment
    infile.close()

#filenames = ["omni15_paired_copynumber_baselineAdjusted_37_RB2012_2015.txt"]
filenames = []
for (_, _, f) in walk(CN_raw_dir):
    filenames += f
    break

outnames = open("ids" + tag + ".txt", "w")
for file in filenames:
    #change 'piv_GC' to 'piv.txt' for non-GC-corrected data.
    if (file.find("omni") == -1 and file.find("Paired") == -1) or file.find("sample_omni") != -1:
        continue
#    if (file.find("REI") == -1):
#        continue
    filebits = file.rstrip().split("_")
    whichomni = filebits[0]
    SNPfile = open(CN_raw_dir + file,"r")
    SNPnames = SNPfile.readline()
    SNPnames = SNPnames.rstrip().split("\t")
    for line in SNPfile:
        line = line.replace("\b", "\t")
        SNPline = line.rstrip().split("\t")
        id = SNPline[0]
        #print id
        idbits = id.rstrip().split("_")
        if (idbits[0].find("-") != -1):
            idbits = idbits[0].split("-")
        if (len(idbits) < 3):
            print("idbits not long enough:", id)
            print("Next entry:", SNPline[1])
            assert(False)
            continue
        id = idbits[0]
        if somepatientsonly and id not in somepatients:
            continue
        sample = idbits[1]
#        if sample != "24007":
#            continue
        if use_canonical:
            if (id + "_" + sample not in whichdata):
                continue
            refomni = whichdata[id + "_" + sample]
            if (refomni != whichomni):
                print("The file " + whichomni + " contains data for " + id + "_" + sample + ", but that data should be found in " + refomni + " instead.")
                continue
        if (len(SNPnames) != len(SNPline)):
            print(str(len(SNPnames)) + " SNPnames and " + str(len(SNPline)) + " SNPs: unequal numbers in file " + file + ", for id '" + id + "_" + sample + "': skipping.")
            for entry in SNPline:
                if (entry.find(":") == -1):
                    x = float(entry)
                    if len(entry) > 13:
                        print(entry)
                        foo()
            continue
#        if (id != "1060"):
#            continue
#        if (sample != "20572"):
#            continue
        outnames.write(SNPline[0] + "\n")
        if path.isfile(CN_out_dir + id + "_" + sample + "_copynumber_all.txt"):
            print("Skipping patient", id, "sample", sample, ": file already exists")
            continue
        print("Writing data for patient", id, ", sample", sample)
        if not checkonly:
            outfile = open(CN_out_dir + id + "_" + sample + "_copynumber_all.txt", "w")
            outfile.write("SNPid\tchr\tpos\tlog2R\n")
            for entry in range(1,len(SNPnames)):
                if SNPnames[entry].find("cnvi") != -1:
                    continue
                label = labels.get(SNPnames[entry])
                if (label == None):
                    #print "No label for", SNPnames[entry]
                    continue
                if (label[0] == "UNKNOWN"):
                    continue;
                try:
                    l2r = float(SNPline[entry])
                except:
                   if SNPline[entry] != "?":
                        print("Non-float value", SNPline[entry], "for", SNPnames[entry])
                   SNPline[entry] = "?"
                if sample in nan_probes and SNPnames[entry] in nan_probes[sample]:
                    #The value was actually faked
                    #print "faked"
                    SNPline[entry] = "?"
                #Unknown value; probably '?'
                outfile.write(SNPnames[entry] + "\t" + label[0] + "\t" + label[1] + "\t" + SNPline[entry] + "\n")
            outfile.close()
    SNPfile.close()
outnames.close()
