# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 13:57:27 2016

@author: lpsmith
"""

from sets import Set
import lucianSNPLibrary as lsl

LOHfilename = "CN_Xiaohong_segmentation/LOH_for_Lucian.txt"
LOHfile = open(LOHfilename, "r")
segs = {}
for line in LOHfile:
    (id, chr, start, end, __, __) = line.split()
    idvec = id.split("_")
    if len(idvec) == 1:
        idvec = id.split("-")
    patient = idvec[0]
    sample = idvec[1]
    if not(patient in segs):
        segs[patient] = {}
    if not sample in segs[patient]:
        segs[patient][sample] = {}
    if not(chr in segs[patient][sample]):
        segs[patient][sample][chr] = []
    segs[patient][sample][chr].append([int(start), int(end), "2\t0"])


#Then read the file Xiaohong sent with the log2R data in it
CNfilename = "CN_Xiaohong_segmentation/copy_number_log2_Lucian.txt"
CNfile = open(CNfilename,"r")

ids = Set()
for line in CNfile:
    (id, chr, start, end, __, ln2r, call) = line.rstrip().split()
    idvec = id.split("_")
    if len(idvec) == 1:
        idvec = id.split("-")
    patient = idvec[0]
    sample = idvec[1]
    if not(patient in segs):
        segs[patient] = {}
    if not sample in segs[patient]:
        segs[patient][sample] = {}
    if not(chr in segs[patient][sample]):
        segs[patient][sample][chr] = []
    if call == "Loss":
        call = "1\t0"
    elif call == "Gain":
        call = "2\t1"
    elif call == "Balanced_gain":
        call = "2\t2"
    elif call == "Double_d":
        call = "0\t0"
    else:
        print "Unknown call", call
    segs[patient][sample][chr].append([int(start), int(end), call])

LOH25Mfile = open("Xiaohong_WGS_segmentation/final_LOH_with_sample_name1.txt", "r")
for line in LOH25Mfile:
    (__, id, __, chr, start, end) = line.rstrip().split()
    idvec = id.split("_")
    if len(idvec) == 1:
        idvec = id.split("-")
    patient = idvec[0]
    sample = idvec[1]
    if not(patient in segs):
        segs[patient] = {}
    if not sample in segs[patient]:
        segs[patient][sample] = {}
    if not(chr in segs[patient][sample]):
        segs[patient][sample][chr] = []
    segs[patient][sample][chr].append([int(start), int(end), "2\t0"])


CNfilename = "Xiaohong_WGS_segmentation/final_copynum_with_sample_name1_NEW1.txt"
CNfile = open(CNfilename,"r")

ids = Set()
for line in CNfile:
    (__, id, __, chr, start, end, __, __, __, call) = line.rstrip().split()
    idvec = id.split("_")
    if len(idvec) == 1:
        idvec = id.split("-")
    patient = idvec[0]
    sample = idvec[1]
    if not(patient in segs):
        segs[patient] = {}
    if not sample in segs[patient]:
        segs[patient][sample] = {}
    if not(chr in segs[patient][sample]):
        segs[patient][sample][chr] = []
    if call == "22" or call=="23":
        call = "1\t0"
    elif call == "32" or "33":
        call = "2\t1"
    elif call == "34":
        call = "2\t2"
    elif call == "41":
        call = "0\t0"
    else:
        print "Unknown call", call
    segs[patient][sample][chr].append([int(start), int(end), call])

for patient in segs:
    for sample in segs[patient]:
        for chr in segs[patient][sample]:
            outvec = segs[patient][sample][chr]
            #outvec = lsl.resegmentLOHes(outvec)
            outfile = open("Xiaohong_segvalidate/" + id + ".txt", "w")
            outfile.write("chr\tstartpos\tendpos\tnA\tnB\n")
            for seg in outvec:
                outfile.write(patient + "\t" + sample + "\t" + str(chr) + "\t" + str(seg[0]) + "\t" + str(seg[1]) + "\t" + seg[2]  + "\n")

#print ids
#print "total number of IDs = " + str(len(ids))
