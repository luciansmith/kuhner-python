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

onlysomepatients = True
somepatients = ["160"]

mutation_file = "snv_plus_indels.20180919.csv"
VAFcompare_file = "VAFcompare.tsv"
outdir = "VAFcompare/"

if not path.isdir(outdir):
    mkdir(outdir)

VAFcompare = open(outdir + VAFcompare_file, "w")
VAFcompare.write("Patient")
VAFcompare.write("\tSampleSet")
VAFcompare.write("\tSampleSet length")
VAFcompare.write("\tSample")
VAFcompare.write("\tSubset length")
VAFcompare.write("\tSubset VAF avg")
VAFcompare.write("\tSubset VAF stdev")
VAFcompare.write("\tFull VAF avg")
VAFcompare.write("\tFull VAF stdev")
VAFcompare.write("\n")


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

def compareVAFs(patient, sampleset, sample, allVAFs, sampleVAFs):
    if len(sampleVAFs)<100:
        return
    VAFcompare.write(patient)
    VAFcompare.write("\t" + str(sampleset))
    VAFcompare.write("\t" + str(len(sampleset)))
    VAFcompare.write("\t" + sample)
    VAFcompare.write("\t" + str(len(sampleVAFs)))
    VAFcompare.write("\t" + str(numpy.average(sampleVAFs)))
    VAFcompare.write("\t" + str(numpy.std(sampleVAFs)))
    VAFcompare.write("\t" + str(numpy.average(allVAFs)))
    VAFcompare.write("\t" + str(numpy.std(allVAFs)))
    VAFcompare.write("\n" )
    

#def writeVAFs(mutations, samplePatientMap):
#    for patient in mutations:
#        for sample in samplePatientMap[patient]:
#            
#        patientVAFs = open(outdir + patient + "_VAFs.tsv", "w")
#        patientVAFs.write("Patient")
#        patientVAFs.write("\tchr")
#        patientVAFs.write("\tpos")
#        for n in range(4):
#            for sample in samplePatientMap[patient]:
#                patientVAFs.write("\t" + sample + "_" + str(n+1) + "_VAF")
#        patientVAFs.write("\tnShared")
#        patientVAFs.write("\n")
#        for chr in mutations[patient]:
#            for pos in mutations[patient][chr]:
#                for alt in mutations[patient][chr][pos]:
#                    patientVAFs.write(patient)
#                    patientVAFs.write("\t" + chr)
#                    patientVAFs.write("\t" + str(pos))
#                    for n in range(4):
#                        for sample in samplePatientMap[patient]:
#                            if len(mutations[patient][chr][pos][alt]) == n+1:
#                                wroteVAF = False
#                                for (msamp, mVAF) in mutations[patient][chr][pos][alt]:
#                                    if msamp==sample:
#                                        patientVAFs.write("\t" + str(mVAF))
#                                        wroteVAF = True
#                                if not wroteVAF:
#                                    patientVAFs.write("\t")
#                            else:
#                                patientVAFs.write("\t")
#                    patientVAFs.write("\n")
#        patientVAFs.close()

def writeVAFs(mutations, samplePatientMap):
    for patient in mutations:
        for sample in samplePatientMap[patient]:
            patientVAFs = open(outdir + patient + "_" + sample + "_VAFs.tsv", "w")
            patientVAFs.write("Patient")
            patientVAFs.write("\tchr")
            patientVAFs.write("\tpos")
            for n in range(6):
                patientVAFs.write("\t" + sample + "_" + str(n+1) + "_VAF")
            patientVAFs.write("\n")
            for chr in mutations[patient]:
                posvec = list(mutations[patient][chr])
                posvec.sort()
                for pos in posvec:
                    for alt in mutations[patient][chr][pos]:
                        outstr = patient + "_" + sample
                        outstr += ("\t" + chr)
                        outstr += ("\t" + str(pos))
                        for n in range(6):
                            if len(mutations[patient][chr][pos][alt]) == n+1:
                                wroteVAF = False
                                for (msamp, mVAF) in mutations[patient][chr][pos][alt]:
                                    if msamp==sample:
                                        outstr += ("\t" + str(mVAF))
                                        wroteVAF = True
                            else:
                                outstr += ("\t")
                        if wroteVAF:
                            patientVAFs.write(outstr + "\n")
            patientVAFs.close()

def writeAllSampleVAFs(mutations, samplePatientMap):
    for patient in mutations:
        patientVAFs = open(outdir + patient + "_allsamples_VAFs.tsv", "w")
        patientVAFs.write("Patient")
        patientVAFs.write("\tchr")
        patientVAFs.write("\tpos")
        patientVAFs.write("\talt")
        for sample in samplePatientMap[patient]:
            patientVAFs.write("\t" + sample)
        for sample in samplePatientMap[patient]:
            patientVAFs.write("\t" + sample + "_present?")
        patientVAFs.write("\n")
        for chr in mutations[patient]:
            posvec = list(mutations[patient][chr])
            posvec.sort()
            for pos in posvec:
                for alt in mutations[patient][chr][pos]:
                    outstr = patient
                    outstr += ("\t" + chr)
                    outstr += ("\t" + str(pos))
                    outstr += ("\t" + alt)
                    wroteVAF = False
                    for sample in samplePatientMap[patient]:
                        for (msamp, mVAF) in mutations[patient][chr][pos][alt]:
                            if msamp==sample:
                                outstr += ("\t" + str(mVAF))
                                wroteVAF = True
                        if not wroteVAF:
                                outstr += ("\t0.0")
                        wroteVAF = False
                    for sample in samplePatientMap[patient]:
                        for (msamp, mVAF) in mutations[patient][chr][pos][alt]:
                            if msamp==sample:
                                outstr += ("\tTrue")
                                wroteVAF = True
                        if not wroteVAF:
                                outstr += ("\tFalse")
                        wroteVAF = False
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
            mutations[patient][chr][pos][alt] = []
        mutations[patient][chr][pos][alt].append((sample, VAF))

writeVAFs(mutations, samplePatientMap)
#writeAllSampleVAFs(mutations, samplePatientMap)

for patient in mutations:
    groups = {}
    all_samples_group = {}
    for chr in mutations[patient]:
        for pos in mutations[patient][chr]:
            for alt in mutations[patient][chr][pos]:
                sampleset = []
                for (sample, VAF) in mutations[patient][chr][pos][alt]:
                    sampleset.append(sample)
                sampleset.sort()
                samp_tuple = tuple(sampleset)
                if samp_tuple not in groups:
                    groups[samp_tuple] = {}
                for (sample, VAF) in mutations[patient][chr][pos][alt]:
                    if sample not in groups[samp_tuple]:
                        groups[samp_tuple][sample] = []
                    groups[samp_tuple][sample].append(VAF)
                    if sample not in all_samples_group:
                        all_samples_group[sample] = []
                    all_samples_group[sample].append(VAF)
    for sampleset in groups:
        for sample in sampleset:
            compareVAFs(patient, sampleset, sample, all_samples_group[sample], groups[sampleset][sample])

VAFcompare.close()
