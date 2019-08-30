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

import lucianSNPLibrary as lsl

onlysomepatients = False
somepatients = ["160"]

mutation_file = "snv_plus_indels.20180919.csv"
outdir = "VAFclusters/"

if not path.isdir(outdir):
    mkdir(outdir)

def getPatientSampleMap():
    patientSampleMap = {}
    samplePatientMap = {}
    callfile = open("calling_evidence.tsv", "r")
    for line in callfile:
        if "Patient" in line:
            continue
        (patient, sample) = line.rstrip().split()[0:2]
        patientSampleMap[sample] = patient
        if patient not in samplePatientMap:
            samplePatientMap[patient] = []
        samplePatientMap[patient].append(sample)
    return patientSampleMap, samplePatientMap

def writeAllSampleVAFs(mutations, samplePatientMap):
    for patient in mutations:
        #Collect a set of clusters
        clustercount = {}
        for chr in mutations[patient]:
            for pos in mutations[patient][chr]:
                for alt in mutations[patient][chr][pos]:
                    cluster = list(mutations[patient][chr][pos][alt].keys())
                    cluster.sort()
                    cluster = tuple(cluster)
                    if cluster not in clustercount:
                        clustercount[cluster] = 0
                    clustercount[cluster] += 1
        clusterlist = []
        for cluster in clustercount:
            if clustercount[cluster] > 40:
                clusterlist.append(cluster)
        for sample in samplePatientMap[patient]:
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


mutations = {}
(patientSampleMap, samplePatientMap) = getPatientSampleMap()
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
        refcnt = int(lvec[-2])
        bafcnt = int(lvec[-1])
        VAF = bafcnt/(refcnt+bafcnt)
        patient = patientSampleMap[sample]
        if onlysomepatients and patient not in somepatients:
            continue
        if patient not in mutations:
            mutations[patient] = {}
        if chr not in mutations[patient]:
            mutations[patient][chr] = {}
        pos = int(pos)
        if pos not in mutations[patient][chr]:
            mutations[patient][chr][pos] = {}
        if alt not in mutations[patient][chr][pos]:
            mutations[patient][chr][pos][alt] = {}
        mutations[patient][chr][pos][alt][sample] = VAF

#writeVAFs(mutations, samplePatientMap)
writeAllSampleVAFs(mutations, samplePatientMap)

