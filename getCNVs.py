#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 13:27:44 2020

@author: lpsmith
"""


from os.path import isfile

def getPatientSampleMap(challenge=False, dipvtet_file=""):
    s2p = {}
    p2s = {}
    if dipvtet_file=="":
        dipvtet_file = "calling_evidence_challenge_inc_odds.tsv"
        if not challenge:
            dipvtet_file = "calling_evidence_odds.tsv"
        
    callfile = open(dipvtet_file, "r")
    for line in callfile:
        if "Patient" in line:
            continue
        lvec = line.rstrip().split()
        (patient, sample) = lvec[0:2]
        ploidy = lvec[-1]
        if ploidy=="Unknown":
            odds = float(lvec[-2])
            if odds > 0.5:
                ploidy = "Diploid"
            else:
                ploidy = "Tetraploid"
        if ploidy=="Diploid":
            ploidy = "diploid"
        if ploidy=="Tetraploid":
            ploidy="tetraploid"
        s2p[sample] = [patient, ploidy]
        if patient not in p2s:
            p2s[patient] = []
        p2s[patient].append(sample)
    return p2s, s2p



def loadDeletionsAndCNVs(samplePatientMap, CNdir="noninteger_processed_CNs/"):
    deletions = {}
    CNVs = {}
    for sample in samplePatientMap:
        (patient, ploidy) = samplePatientMap[sample]
        filename = CNdir + patient + "_" + sample + "_g500_" + ploidy + "_nonint_CNs.txt"
        if not isfile(filename):
            filename = CNdir + patient + "_" + sample + "_g550_" + ploidy + "_nonint_CNs.txt"
        for line in open(filename, "r"):
            if "patient" in line:
                continue
            lvec = line.rstrip().split()
            if lvec[7] == "0" or lvec[8] == "0":
                (chrom, start, end) = lvec[2:5]
                start = int(start)
                end = int(end)
                if patient not in deletions:
                    deletions[patient] = {}
                if sample not in deletions[patient]:
                    deletions[patient][sample] = {}
                if chrom not in deletions[patient][sample]:
                    deletions[patient][sample][chrom] = []
                deletions[patient][sample][chrom].append((start, end))
            intA = lvec[7]
            intB = lvec[8]
            if intA=="NA" or intB=="NA":
                continue
            intA = int(intA)
            intB = int(intB)
            if intB<intA:
                temp = intA
                intA = intB
                intB = temp
            call = (intA, intB)
            (chrom, start, end) = lvec[2:5]
            start = int(start)
            end = int(end)
            if patient not in CNVs:
                CNVs[patient] = {}
            if sample not in CNVs[patient]:
                CNVs[patient][sample] = {}
            if chrom not in CNVs[patient][sample]:
                CNVs[patient][sample][chrom] = []
            CNVs[patient][sample][chrom].append((start, end, call))
    return deletions, CNVs


def getCNVCall(patient, sample, chrom, pos, CNVs):
    if patient not in CNVs:
        assert(False)
        return (-1, -1)
    if sample not in CNVs[patient]:
        assert(False)
        return (-1, -1)
    if chrom not in CNVs[patient][sample]:
        print("No chromosome", str(chrom), "found.")
        assert(False)
        return (-1, -1)
    for (start, end, call) in CNVs[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return call
    return (-1, -1)

patientSampleMap, samplePatientMap = getPatientSampleMap()
deletions, CNVs = loadDeletionsAndCNVs(samplePatientMap)
print(str(getCNVCall("995", "24945", "9", 500000, CNVs)))