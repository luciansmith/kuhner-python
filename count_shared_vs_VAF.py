#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:42:05 2019

@author: lpsmith
"""


import csv

import lucianSNPLibrary as lsl

mutation_file = "snv_plus_indels.twoPlus.20181030.csv"


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


mutations = {}
with open(mutation_file, 'r') as csvfile:
    for lvec in csv.reader(csvfile):
        if "DNANum" in lvec[0]:
            continue
        (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
        if (is_2p=="f"):
            continue
        ref = int(lvec[-2])
        alt = int(lvec[-1])
        VAF = alt/(ref+alt)
        
    #    if ("N" in sample):
    #        continue
        if sample not in mutations:
            mutations[sample] = {}
        if chr not in mutations[sample]:
            mutations[sample][chr] = {}
        pos = int(pos)
        mutations[sample][chr][pos] = (ref, alt, VAF)

(patientSampleMap, samplePatientMap) = getPatientSampleMap()

shared_high = 0
shared_low = 0
shared_both = 0
private_high = 0
private_low = 0

for patient in samplePatientMap:
    patientmuts = {}
    for sample in samplePatientMap[patient]:
        for chr in mutations[sample]:
            if chr not in patientmuts:
                patientmuts[chr] = {}
            for pos in mutations[sample][chr]:
                if pos not in patientmuts[chr]:
                    patientmuts[chr][pos] = {}
                patientmuts[chr][pos][sample] = mutations[sample][chr][pos]
    for chr in patientmuts:
        for pos in patientmuts[chr]:
            samples = list(patientmuts[chr][pos].keys())
            canon_alt = patientmuts[chr][pos][samples[0]][1]
            nhigh = 0
            nlow = 0
            for sample in samples:
                if patientmuts[chr][pos][sample][1] == canon_alt:
                    if patientmuts[chr][pos][sample][2] >= 0.25:
                        nhigh += 1
                    else:
                        nlow += 1
            if nhigh + nlow == 1:
                if nhigh==1:
                    private_high += 1
                else:
                    private_low += 1
            else:
                if nhigh > 0 and nlow > 0:
                    shared_both +=1
                elif nhigh > 0:
                    shared_high += 1
                else:
                    shared_low += 1

print(str(shared_high), "shared high")
print(str(shared_low), "shared low")
print(str(shared_both), "shared both")
print(str(private_high), "private high")
print(str(private_low), "private low")
