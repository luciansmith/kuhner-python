# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:43:41 2018

@author: Lucian
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy
from ete3 import Tree

import numpy
import math
import matplotlib.pyplot as plt
import csv

# replace with call to imp.source_load, renamed "lsl" as "lps"
#import lucianSNPLibrary as lsl

onlysomepatients = False
somepatients = ["160"]

#mutation_file = "snv_plus_indels.twoPlus.20181030.csv"
mutation_file = "lucian_from_kanika.csv"
outdir = "VAFclusters/"

if not path.isdir(outdir):
    mkdir(outdir)

def writeAllSampleVAFs(mutations, patientSampleMap, deletions):
    for patient in mutations:
        #Collect a set of clusters
        print("Writing data for patient", patient)
        clustercount = {}
        for chr in mutations[patient]:
            for pos in mutations[patient][chr]:
                for alt in mutations[patient][chr][pos]:
                    cluster = list(mutations[patient][chr][pos][alt].keys())
                    cluster.sort()
                    for sample in patientSampleMap[patient]:
                        if sample in cluster:
                            continue
                        if isDeleted(patient, sample, chr, pos, deletions):
                            cluster.append("-" + sample)
                            mutations[patient][chr][pos][alt]["-" + sample] = ""
                    cluster = tuple(cluster)
                    if cluster not in clustercount:
                        clustercount[cluster] = 0
                    clustercount[cluster] += 1
        clusterlist = []
        for cluster in clustercount:
#            if clustercount[cluster] > 40:
                clusterlist.append(cluster)
        for sample in patientSampleMap[patient]:
            patientVAFs = open(outdir + patient + "_" + sample + "_VAFs.tsv", "w")
            patientVAFs.write("Patient")
            patientVAFs.write("\tSample")
            patientVAFs.write("\tchr")
            patientVAFs.write("\tpos")
            patientVAFs.write("\talt")
            for cluster in clusterlist:
                if sample in cluster:
                    patientVAFs.write("\t" + str(cluster))
            patientVAFs.write("\n")
            
            for chr in mutations[patient]:
                posvec = list(mutations[patient][chr].keys())
                posvec.sort()
                for pos in posvec:
                    for alt in mutations[patient][chr][pos]:
                        outstr = patient
                        outstr += ("\t" + sample)
                        outstr += ("\t" + chr)
                        outstr += ("\t" + str(pos))
                        outstr += ("\t" + alt)
                        writeVAF = False
                        for cluster in clusterlist:
                            if sample in cluster and sample in mutations[patient][chr][pos][alt]:
                                #We have to write something: VAF or a tab.
                                outstr += "\t"
                                if set(cluster) == set(mutations[patient][chr][pos][alt].keys()):
                                    outstr += str(mutations[patient][chr][pos][alt][sample])
                                    writeVAF = True
                        if writeVAF:
                            patientVAFs.write(outstr + "\n")
        patientVAFs.close()

def isDeleted(patient, sample, chrom, pos, deletions):
    if patient not in deletions:
        return False
    if sample not in deletions[patient]:
        return False
    if chrom not in deletions[patient][sample]:
        return False
    for (start, end) in deletions[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return True
    return False

def getSampleStatuses():
    statuses = {}
    with open("P01CA91955-WGS80-Full-Pilot-Samples.csv", "r") as csvfile:
        for lvec in csv.reader(csvfile):
            if "RandomID" in lvec:
                continue
            sample = lvec[6]
            statuses[sample] = {}
            statuses[sample]["patient"] = lvec[0]
            statuses[sample]["prog"] = lvec[1]
            statuses[sample]["gender"] = lvec[2]
            time = lvec[8]
            if time=="LongFollowUp":
                time = "T3"
            if time=="Index":
                time = "T2"
            if sample == "23521":
                time = "T2"
            statuses[sample]["time"] = time
    return statuses


import imp
#import lucianSNPLibrary as lps
lps = imp.load_source("lps","/home/mkkuhner/Papers/phylo/lucianSNPLibrary.py")
mutations = {}
(__, samplePatientMap) = lps.getPatientSampleMap()

#Delete all samples except T1 and T2.
statuses = getSampleStatuses()
for sample in statuses:
    if statuses[sample]["time"] not in ["T1", "T2"]:
        if sample in samplePatientMap:
            del samplePatientMap[sample]
        
        

patientSampleMap = {}

with open(mutation_file, 'r') as csvfile:
    for lvec in csv.reader(csvfile):
        if "DNANum" in lvec[0]:
            continue
        (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
        if (is_snv=="f"):
            continue
        if (is_2p=="f"):
            continue
    #    if ("N" in sample):
    #        continue
        if sample not in samplePatientMap:
            continue
        refcnt = int(lvec[-2])
        bafcnt = int(lvec[-1])
        VAF = bafcnt/(refcnt+bafcnt)
        patient = samplePatientMap[sample][0]
        if patient not in patientSampleMap:
            patientSampleMap[patient] = set()
        patientSampleMap[patient].add(sample)
        if onlysomepatients and patient not in somepatients:
            continue
        if patient not in mutations:
            mutations[patient] = {}
            print("Reached patient", patient)
        if chr not in mutations[patient]:
            mutations[patient][chr] = {}
        pos = int(pos)
        if pos not in mutations[patient][chr]:
            mutations[patient][chr][pos] = {}
        if alt not in mutations[patient][chr][pos]:
            mutations[patient][chr][pos][alt] = {}
        mutations[patient][chr][pos][alt][sample] = VAF

deletions = lps.loadDeletions(samplePatientMap)
writeAllSampleVAFs(mutations, patientSampleMap, deletions)

