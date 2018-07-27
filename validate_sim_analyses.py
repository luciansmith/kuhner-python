#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

simdir = "siminput/"

balanced_outdir = "simoutput/"

bafrawdata = {}
patient_samples = {}

onlysomepatients = False
somepatients = ["572"]
onlysomechroms = False
somechroms = ["9"]

if not path.isdir(balanced_outdir):
    mkdir(balanced_outdir)

ignore_cnvis = True

bafWtLow = 0.4
bafWtHigh = 0.65

#original cutoff: 0.35 to 0.65, or 0.15 out from 0.5.  The above values account for dye bias.

def readBafNormal(patient):
    bafrawdata = {}
    bafwt = {}
#    allcnvis = []
    bafnormal = readlink(simdir + patient + "_Normal_BAF.txt")
    if not(isfile(bafnormal)):
        print("ERROR:  no Normal BAF file found for patient", patient)
        return ({}, {})
    bafnormal = open(bafnormal, "r")
    print("Reading BAF normal data for patient", patient)
    for line in bafnormal:
        lvec = line.split()
        if line.find("Chr") != -1:
            continue
        if ignore_cnvis and line.find("cnvi") != -1:
            continue
        try:
            value = float(lvec[3])
        except:
            continue
#        if (line.find("cnvi") != -1):
#            allcnvis.append(value)
        if (value < bafWtLow or value > bafWtHigh):
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if onlysomechroms and chr not in somechroms:
            continue
        if chr not in bafrawdata:
            bafrawdata[chr] = {}
        if chr not in bafwt:
            bafwt[chr] = {}
        bafrawdata[chr][pos] = {}
        bafwt[chr][pos] = value
    bafnormal.close()
#    print("All normal CNVI bafs:")
#    lsl.createPrintAndSaveHistogram(allcnvis, "", .01)
    return bafrawdata, bafwt

def readBafSamples(patient, sample, bafrawdata):
    labels = []
    all_samples = []
    baffile = readlink(simdir + patient + "_" + sample + "_BAF.txt")
    if not(isfile(baffile)):
        print("ERROR:  no BAF file found for patient", patient, "sample", sample)
        bafrawdata = {}
        return
    print("Reading BAF sample data for patient", patient, "sample", sample)
    baffile = open(baffile, "r")
    for line in baffile:
        lvec = line.split()
        if "Chr" in line:
            labels = lvec
            for l in range(2,len(labels)):
                all_samples.append(labels[l].split('"')[1])
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if onlysomechroms and chr not in somechroms:
            continue
        if chr not in bafrawdata:
            continue
        if pos not in bafrawdata[chr]:
            continue
        for p in range(3,len(lvec)):
            sample = labels[p-1].split('"')[1]
            try:
                bafrawdata[chr][pos][sample] = float(lvec[p])
                if line.find("cnvi") != -1:
                    allbafs["cnvi"].append(float(lvec[p]))
                else:
                    bafrawdata[chr][pos][sample] = float(lvec[p])
                    allbafs["normal"].append(float(lvec[p]))
            except:
                continue
    baffile.close()
    return allbafs, all_samples

def getIsegsFromCopynumberFileFor(patient):
    #This function reads the lowest-gamma-value copynumber file it can, reads it, and returns it as 'isegs'.
    isegs = {}
    #glist = ("100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000")
    glist = list("1")
    glist[0] = "500"
    gindex = 0
    gamma = glist[gindex]

    root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"
    isegfilename = root_dir + patient + "_copynumber_segments.txt"
    while not(isfile(isegfilename)) and gindex<len(glist):
        gamma = glist[gindex]
        root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"
        isegfilename = root_dir + patient + "_copynumber_segments.txt"
        gindex += 1
    if gindex > len(glist):
        print("Cannot find any copynumber file for patient", patient)
        return {}
    isegfile = open(isegfilename, "r")
    for line in isegfile:
        if line.find("Chr") != -1:
            continue
        (chr, start, end, nlogr, nbaf) = line.split()
        if onlysomechroms and chr not in somechroms:
            continue
        nbaf = int(nbaf)
#        if nbaf <10:
#            continue
        start = int(start)
        end = int(end)
        if chr=="3" and end==198837449:
            end = 198022430
        elif chr=="19" and end==121485079:
            end = 59128983
        elif chr=="21" and (end==106115710 or end==207101487):
            end = 48129895
        elif chr=="22" and end==51666786:
            end = 51304566
        if chr not in isegs:
            isegs[chr] = []
        isegs[chr].append([start, end, [], [], {}])
    return isegs

def readPloidyFile(ploidyfile):
    ploidies = {}
    pf = open(ploidyfile)
    for line in pf:
        if line.find('x') != -1:
            continue
        (id, val) = line.split()
        sample = id.split('"')[1].split('_')[1]
        ploidies[sample] = val
    pf.close()
    return ploidies

def readPurityFile(purityfile):
    purities = {}
    pf = open(purityfile)
    for line in pf:
        if line.find('x') != -1:
            continue
        (id, val) = line.split()
        sample = id.split('"')[1].split('_')[1]
        purities[sample] = val
    pf.close()
    return purities

def readSegmentationFile(ascsegfile, all_samples):
    totsca = {}
    totsca["overall"] = set()
    for sample in all_samples:
        totsca[sample] = set()
    osegs = {}
    for line in open(ascsegfile):
        if line.find("Chr") != -1:
            continue
        (__, sample, chr, start, end, nprobes, mbaf, logr, nA, nB) = line.split()
        if onlysomechroms and chr not in somechroms:
            continue
        sample = sample.split('"')[1]
        start = int(start)
        end = int(end)
        if (nA == nB):
            continue
        if sample not in osegs:
            osegs[sample] = {}
        if chr not in osegs[sample]:
            osegs[sample][chr] = []
        osegs[sample][chr].append((start, end))
        totsca[sample].add((chr, start, end))
        totsca["overall"].add((chr, start, end))
    return (osegs, totsca)

def addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca):
    if onlysomechroms and chr not in somechroms:
        return
    if chr == "":
        return
    if chr == "23":
        return
    if chr == "24":
        return
    if chr=="0":
        #print("There's a chromosome zero segment, weirdly:", full_sample, chr, start, end)
        return
    svec = full_sample.split("-")
    if len(svec) < 4:
        svec = full_sample.split("_")
    patient = svec[0]
    if onlysomepatients and patient not in somepatients:
        return
    sample = svec[0] + "_" + svec[1]
    if patient not in Xiaohong_segments:
        Xiaohong_segments[patient] = {}
        totsca[patient] = {}
    if sample not in Xiaohong_segments[patient]:
        Xiaohong_segments[patient][sample] = {}
        totsca[patient][sample] = set()
    if chr not in Xiaohong_segments[patient][sample]:
        Xiaohong_segments[patient][sample][chr] = []
    start = int(start)
    end = int(end)
    Xiaohong_segments[patient][sample][chr].append((start, end))
    totsca[patient][sample].add((chr, start, end))
    if "overall" not in totsca[patient]:
        totsca[patient]["overall"] = set()
    totsca[patient]["overall"].add((chr, start, end))
    #print("Adding", chr, start, end)

def readXiaohongWGSLOHFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (__, full_sample, __, chr, start, end) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca)
    xfile.close()

def readXiaohong1MLOHFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    lastseg = ["", 0, 0, ""]
    for line in xfile:
        (full_sample, chr, start, end, __, __) = line.split()
        start =int(start)
        end = int(end)
        if full_sample == lastseg[3] and chr == lastseg[0] and start-lastseg[2] < 5000:
            print("Combining", chr, str(lastseg[1]), str(lastseg[2]), "with", start, end)
            lastseg[2] = end
        else:
            addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)
            lastseg = [chr, start, end, full_sample]
    xfile.close()
    addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)

def readXiaohongCopynumFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    lastseg = ["", 0, 0, 0, ""]
    for line in xfile:
        lvec = line.split()
        if len(lvec) == 7:
            (full_sample, chr, start, end, __, __, code) = line.split()
        elif len(lvec) == 10:
            (__, full_sample, __, chr, start, end, __, __, __, code) = line.split()
        else:
            print("Incorrect line length:", line)
            continue
        start =int(start)
        end = int(end)
        if code == "Double_d" or code=="Balanced_gain" or code == "34" or code=="41":
            #balanced gain or loss: continue
            addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)
            lastseg = ["", 0, 0, 0, ""]
            continue
        if full_sample == lastseg[4] and chr == lastseg[0] and code == lastseg[3] and start-lastseg[2] < 5000:
            #print("Combining", chr, str(lastseg[1]), str(lastseg[2]), str(lastseg[3]), "with", start, end, code)
            lastseg[2] = end
        else:
            addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)
            lastseg = [chr, start, end, code, full_sample]
    xfile.close()
    addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)

def readAllXiaohongSegmentation():
    print("Reading all Xiaohong segmentation.")
    Xiaohong_segments = {}
    totsca = {}
    files = []
    for (__, __, f) in walk(Xdir_WGS):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("test") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohongWGSLOHFile(Xdir_WGS + f, Xiaohong_segments, totsca)
        else:
            readXiaohongCopynumFile(Xdir_WGS + f, Xiaohong_segments, totsca)
    files = []
    for (__, __, f) in walk(Xdir_1M):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohong1MLOHFile(Xdir_1M + f, Xiaohong_segments, totsca)
        else:
            readXiaohongCopynumFile(Xdir_1M + f, Xiaohong_segments, totsca)
    return Xiaohong_segments, totsca

def storeMatchesInIsegs(isegs, bafrawdata):
    for chr in isegs:
        if onlysomechroms and chr not in somechroms:
            continue
        for iseg in isegs[chr]:
            patterns = {}
            for pos in bafrawdata[chr]:
                if pos >= iseg[0] and pos <= iseg[1]:
                    rdsamples = list(bafrawdata[chr][pos].keys())
                    rdsamples.sort()
                    for sample1 in range(0, len(rdsamples)-1):
                        for sample2 in range(sample1+1, len(rdsamples)):
                            s1 = rdsamples[sample1]
                            s2 = rdsamples[sample2]
                            try:
                                val1 = bafrawdata[chr][pos][s1]
                                val2 = bafrawdata[chr][pos][s2]
                            except:
                                continue
                            segpair = (s1, s2)
                            if segpair not in patterns:
                                patterns[segpair] = [0, 0]
                            if (val1 > 0.5 and val2 > 0.5) or (val1 < 0.5 and val2 < 0.5):
                                patterns[segpair][1] += 1
                            patterns[segpair][0] += 1
            for segpair in patterns:
                pattern = patterns[segpair]
                if pattern[0] < 10:
                    continue
                ratio = pattern[1]/pattern[0]
                if ratio > 0.95 or ratio < 0.05:
                    iseg[2].append(segpair)
                else:
                    iseg[3].append(segpair)

def scoreAnalysis(isegs, osegs, sample, analysis):
    for chr in isegs:
        if onlysomechroms and chr not in somechroms:
            continue
        for iseg in isegs[chr]:
            if sample not in iseg[4]:
                iseg[4][sample] = {}
            istart = iseg[0]
            iend = iseg[1]
            matches = iseg[2]
            nonmatches = iseg[3]
            sample_in_matches = False
            sample_in_nonmatches = False
            for match in matches:
                if sample in match:
                    sample_in_matches = True
                    break
            for nonmatch in nonmatches:
                if sample in nonmatch:
                    sample_in_nonmatches = True
                    break
            call = "Balanced"
            if sample in osegs and chr in osegs[sample]:
                for oseg in osegs[sample][chr]:
                    if oseg[0] <= istart and oseg[1] >= iend:
                        call = "Unbalanced"
                        break
                    elif oseg[0] > istart and oseg[0] < iend:
                        call = "None"
                    elif oseg[1] > istart and oseg[1] < iend:
                        call = "None"
                    #Xiaohong's segments sometimes overcall an area as both something and CNLOH, so 'unbalanced' should override 'None'.
            if call=="None":
                iseg[4][sample][analysis] = "NC" #"No call"
                continue
            if not sample_in_matches and not sample_in_nonmatches:
                #segment might be too short: sample not present
                iseg[4][sample][analysis] = "S" #"Too short?"
                continue
            other_matches = False
            for (m1, m2) in matches:
                for nonmatch in nonmatches:
                    if sample in nonmatch and (m1 in nonmatch or m2 in nonmatch):
                        other_matches = True
                        break
                    
            if call == "Unbalanced":
                if sample_in_matches:
                    iseg[4][sample][analysis] = "TP" #"True positive"
                elif other_matches:
                    iseg[4][sample][analysis] = "FP" #"False positive"
                else:
                    iseg[4][sample][analysis] = "UP" #"Unvalidatable positive"
            elif call == "Balanced":
                if sample_in_matches:
                    iseg[4][sample][analysis] = "FN" #"False negative"
                elif other_matches:
                    iseg[4][sample][analysis] = "TN" #"True negative"
                else:
                    iseg[4][sample][analysis] = "UN" #"Unvalidatable negative"

def writeBalanceFiles(isegs, patient, all_samples):
    for sample in all_samples:
        outfile = open(balanced_outdir + sample + "_balanced_calls.tsv", "w")
        outfile.write("Patient\tSample\tChr\tStart\tEnd\tCall\n")
        for chr in isegs:
            for iseg in isegs[chr]:
                outfile.write(patient)
                outfile.write("\t" + sample.split("_")[1])
                outfile.write("\t" + chr)
                outfile.write("\t" + str(iseg[0]))
                outfile.write("\t" + str(iseg[1]))
                matches = iseg[2]
                nonmatches = iseg[3]
                sample_in_matches = False
                sample_in_nonmatches = False
                for match in matches:
                    if sample in match:
                        sample_in_matches = True
                        break
                for nonmatch in nonmatches:
                    if sample in nonmatch:
                        sample_in_nonmatches = True
                        break
                other_matches = False
                for (m1, m2) in matches:
                    for nonmatch in nonmatches:
                        if sample in nonmatch and (m1 in nonmatch or m2 in nonmatch):
                            other_matches = True
                            break
                if sample_in_matches:
                    outfile.write("\tUnbalanced")
                elif not sample_in_nonmatches:
                    outfile.write("\tShort")
                elif other_matches:
                    outfile.write("\tBalanced")
                else:
                    outfile.write("\tUnknown")
                outfile.write("\n")

def writeDerivedStatistic(ov_out, ov, sample, analysis, numerators, denoms):
    numerator = 0
    denom = 0
    for num in numerators:
        numerator += ov[sample][analysis][num]
        denom += ov[sample][analysis][num]
    for den in denoms:
        denom += ov[sample][analysis][den]
    value = 0
    if denom != 0:
        value = numerator/denom
        ov_out.write("\t" + str(value))
    else:
        ov_out.write("\t--")
    return value


def writeSummary(isegs, patient, all_samples, all_analyses):
    if not onlysomepatients and isfile(outdir + patient + "_analysis_overview.txt"):
        print("Skipping patient", patient, ": analysis already exists.")
        return
    overview = {}
    overview_bases = {}
    calls = ["TP", "FP", "UP", "TN", "FN", "UN", "NC", "S"]
    for sample in all_samples:
        overview[sample] = {}
        overview_bases[sample] = {}
        for analysis in all_analyses:
            overview[sample][analysis] = {}
            overview_bases[sample][analysis] = {}
            for call in calls:
                overview[sample][analysis][call] = 0
                overview_bases[sample][analysis][call] = 0
    summary_out = open(outdir + patient + "_full_analysis.tsv", "w")
    summary_out.write("Patient")
    summary_out.write("\tSample")
    summary_out.write("\tchr")
    summary_out.write("\tstart")
    summary_out.write("\tend")
    for analysis in all_analyses:
        summary_out.write("\t" + analysis)
    summary_out.write("\n")
    for sample in all_samples:
        for chr in isegs:
            for iseg in isegs[chr]:
                summary_out.write(patient)
                summary_out.write("\t" + sample.split("_")[1])
                summary_out.write("\t" + chr)
                summary_out.write("\t" + str(iseg[0]))
                summary_out.write("\t" + str(iseg[1]))
                for analysis in all_analyses:
                    summary_out.write("\t" + iseg[4][sample][analysis])
                    overview[sample][analysis][iseg[4][sample][analysis]] += 1
                    overview_bases[sample][analysis][iseg[4][sample][analysis]] += iseg[1] - iseg[0]
                summary_out.write("\n")
    summary_out.close()

    overview_out = open(outdir + patient + "_analysis_overview.tsv", "w")
    overview_out.write("Patient")
    overview_out.write("\tSample")
    overview_out.write("\tAnalysis gamma")
    overview_out.write("\tAnalysis ploidy")
    for call in calls:
        overview_out.write("\t" + call + "_num")
    for call in calls:
        overview_out.write("\t" + call + "_len")
    overview_out.write("\tnum_accuracy")
    overview_out.write("\tlen_accuracy")
    overview_out.write("\n")
    for sample in all_samples:
        for analysis in all_analyses:
            overview_out.write(patient)
            overview_out.write("\t" + sample.split("_")[1])
            analysis_split = -1
            if "tetraploid" in analysis:
                analysis_split = analysis.find("t")
                atype = "tetraploid"
            if "diploid" in analysis:
                analysis_split = analysis.find("d")
                atype = "diploid"
            elif "eight" in analysis:
                analysis_split = analysis.find("e")
                atype = "eight"
            if analysis_split == -1:
                overview_out.write("\t0\tXiaohong")
                analysis_split = len(analysis)
                atype = "Xiaohong"
            else:
                overview_out.write("\t" + analysis[:analysis_split] + "\t" + analysis[analysis_split:])
            for call in calls:
                overview_out.write("\t" + str(overview[sample][analysis][call]))
            for call in calls:
                overview_out.write("\t" + str(overview_bases[sample][analysis][call]))
            writeDerivedStatistic(overview_out, overview, sample, analysis,  ["TP", "TN"], ["FP", "FN"])
            writeDerivedStatistic(overview_out, overview_bases, sample, analysis,  ["TP", "TN"], ["FP", "FN"])
            overview_out.write("\n")
    overview_out.close()


(Xiaohong_segments, X_totsca) = readAllXiaohongSegmentation()

files = []
for (__, __, f) in walk(BAF_dir):
    files += f
for f in files:
    if f.find("_Normal_BAF.txt") == -1:
        continue
    patient = f.split("_")[0]
    if (onlysomepatients and patient not in somepatients):
        continue

    bafrawdata, bafwt = readBafNormal(patient)
    if (len(bafrawdata)==0):
        continue
    allbafs, all_samples = readBafSamples(patient, bafrawdata)

    isegs = getIsegsFromCopynumberFileFor(patient)
    storeMatchesInIsegs(isegs, bafrawdata)

    ploidies = {}
    purities = {}
    all_analyses = ["Xiaohong"]

    for gamma in gamma_list:
        print("Processing results from a gamma of", gamma)
        root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"
        ploidies[gamma] = {}
        purities[gamma] = {}

        for subdir in subdirs:
            ploidyfile = root_dir + subdir + "/" + patient + "_fcn_ascat_ploidy.txt"
            if not(isfile(ploidyfile)):
                print("No ploidy for", patient, "gamma", gamma, "for the", subdir, "constrained run: ASCAT failure for this patient.")
                continue
            ploidies[gamma][subdir] = readPloidyFile(ploidyfile)
            purityfile = root_dir + subdir + "/" + patient + "_fcn_ascat_cont.txt"
            if not(isfile(purityfile)):
                print("No purity for", patient, "gamma", gamma, "for the", subdir, "constrained run: ASCAT failure for this patient.")
                continue
            purities[gamma][subdir] = readPurityFile(purityfile)
            ascsegfile = root_dir + subdir + "/" + patient + "_fcn_ascat_segments.txt"
            if not(isfile(ascsegfile)):
                print("ASCAT seems to have failed for", root_dir, ",",subdir,",",patient,".")
                continue

            all_analyses.append(gamma+subdir)
            (osegs, totsca) = readSegmentationFile(ascsegfile, all_samples)
            for sample in all_samples:
                scoreAnalysis(isegs, osegs, sample, gamma+subdir)

    print("Running SCA analysis for Xiaohong data for patient", patient)
    for sample in all_samples:
        scoreAnalysis(isegs, Xiaohong_segments[patient], sample, "Xiaohong")

    writeSummary(isegs, patient, all_samples, all_analyses)
    writeBalanceFiles(isegs, patient, all_samples)


