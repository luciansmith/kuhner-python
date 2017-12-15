# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk

import lucianSNPLibrary as lsl
import string
import numpy


nsamples_min = 10 #Arbitrary value: minimum number of samples we require
#nsamples_max = 12

# read the file that correlates patient data with which omni file the 'canonical' version of that data can be found.
infilename = "CN_raw_data/20160621_sample_omni.txt"
infile = open(infilename,"r")

whichdata = {}
for line in infile:
    if (line[0]=="PatientID"):
        continue
    line = line.rstrip().split("\t")
    (id, num, experiment) = line
    whichdata[id + "_" + num] = experiment
infile.close()

#fullseg_filenames = ["three_formal_cy_omni_15_b37RB.txt"]
fullseg_filenames = []
for (_, _, f) in walk("full_segmentation/"):
    fullseg_filenames += f
    break        

def empty22():
    return [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], ]

patientdata = []
for filename in fullseg_filenames:
    if (string.find(filename, "omni") == -1):
        continue
    filebits = filename.rstrip().split("_") #Format is "three_formal_cy_omni_xx_b37RB.txt" where xx is the omni number.
    omniid = "omni" + filebits[4]
    omniid = omniid.replace("mix", "Mix")
    print "Now processing omni file", omniid
    patientlist = set()
    for patients in whichdata.items():
        if (patients[1] == omniid):
            patientlist.add(patients[0])
    handle = open("full_segmentation/" + filename, "r")
    patientdata = []
    onepatientdata = empty22()
    patients = []
    patient = ""
    lastpatient = ""
    lastchr = 0
    lastline = []
    for line in handle:
        (chr, start, end, id, classification, bps, mean, nmarkers, pval) = line.rstrip().split("\t")
        if (chr == "chromosome"):
            continue
        chr = chr[3:]
        try:
            chr = int(chr)
        except:
            continue
        id = id.split("_")
        patient = id[0] + "_" + id[1]
        if patient not in patientlist:
            continue
        start = int(start)
        end = int(end)
        if patient != lastpatient and lastpatient != "":
            while(lastchr < 22):
                onepatientdata[lastchr].append([0, "inf", "?", "?", []])
                lastchr += 1
            patientdata.append(onepatientdata)
            onepatientdata = empty22()
            lastchr = 0
            lastline = []
            patients.append(lastpatient)

        lastpatient = patient
        if (chr != lastchr):
            if len(lastline) > 0:
                onepatientdata[lastchr].append([lastline[1], "inf", "?", "?", []])
            lastchr += 1
            while(lastchr < chr):
                onepatientdata[lastchr].append([0, "inf", "?", "?", []])
                lastchr += 1
            lastchr = chr #Just to make sure.
            if (start != 0):
                onepatientdata[chr].append([0, start, "?", "?", []])
        elif lastline[1] != start:
            onepatientdata[chr].append([lastline[1], start, "?", "?", []])
                
        lastline = [start, end, mean, nmarkers, []]
        onepatientdata[chr].append(lastline)
        
    while(lastchr < 22):
        onepatientdata[lastchr].append([0, "inf", "?", "?", []])
        lastchr += 1
    patientdata.append(onepatientdata)
    patients.append(lastpatient)
    handle.close()

    for p in range(0, len(patients)):
        patient = patients[p]
        try:
            handle = open("CN_filtered_data/" + patient + "_copynumber_all.txt", "r")
        except:
            print "No such file CN_filtered_data/" + patient + "_copynumber_all.txt"
            continue                
        for line in handle:
            (id, chr, pos, log2r) = line.rstrip().split("\t")
            if (id == "SNPid"):
                continue
            chr = int(chr)
            if (chr==0):
                continue
            if (chr>=23):
                continue
            pos = int(pos)
            if (pos==0):
                continue
            try:
                log2r = float(log2r)
            except:
                continue                    
            for segment in patientdata[p][chr]:
                if (pos <= segment[0]): #The 'start' is *excluded* from the segment.
                    continue
                if (pos >  segment[1]): #The 'end' is *included* in the segment
                    continue
                segment[4].append(log2r)
                break
        handle.close()
        #Now print out the data:
        handle = open("full_segmentation_output/" + patient + "_segmented.txt", "w")
        handle.write("chr\tstart\tend\tpartekMean\tpartekNMarkers\tnmarkers\tmeanlog2r\n")
        for c in range(1,23):
            for segment in patientdata[p][c]:
                (start, end, pMean, pNMarkers, data) = segment
                outline = str(c) + "\t" + str(start) + "\t" + str(end) + "\t" + pMean + "\t" + pNMarkers + "\t" + str(len(data)) + "\t" + str(numpy.mean(data)) + "\n"
                handle.write(outline)
                if pNMarkers != "?" and int(pNMarkers) != len(data):
                    print "Anomaly for", patient, c, start, end, ":  Partek nmarkers = ", pNMarkers, ", len data =", len(data)
        handle.close()
        print "Wrote data for", patient
        
            
                
