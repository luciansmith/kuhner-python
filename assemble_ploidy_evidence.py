#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 14:09:13 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os.path import isfile
from os import mkdir

#Process various input files to create an evidence file

initcall = "initial_calling_evidence.tsv"
finalcall_file = "tom_final_DNA_calls.tsv"
flow_summary = "flow_summary.txt"
call_summary = "Xiaohong_pASCAT_compare/xiaocompare_summary.tsv"
goodness_dir = "gamma_test_output/pASCAT_input_g"

include_challenge = True

if include_challenge:
    outfile = "calling_evidence_challenge_inc.tsv"
else:
    outfile = "calling_evidence.tsv"

def writeNewLine(f, patient, vals):
    f.write(patient)
    for val in vals:
        f.write("\t" + str(val))
    f.write("\n")


def readFlowSummary(flow):
    flowdata = {}
    for line in open(flow, "r"):
        if "Patient" in line:
            continue
        lvec = line.split()
        (progression, patient, ratio) = lvec[0:3]
        aneuploid_strs = lvec[3:len(lvec)]
        aneuploids = set()
        for aneuploid in aneuploid_strs:
            aneuploids.add(float(aneuploid))
        ratio = ratio.replace(':','::')
        flowdata[patient] = (ratio, progression, aneuploids)
    return flowdata

def readCallSummary(call):
    calldata = {}
    for line in open(call, "r"):
        if "Patient" in line:
            continue
        (patient, sample, gamma, ploidy, ploidyval, purval, __, __, __, __, __, __, __, __, __, __, Xacc, Aacc) = line.split()
        ploidyval = float(ploidyval)
        purval = float(purval)
        if Aacc == "--":
            Aacc = 0.0
        else:
            Aacc = float(Aacc)
        if patient not in calldata:
            calldata[patient] = {}
        if sample not in calldata[patient]:
            calldata[patient][sample] = {}
        calldata[patient][sample][ploidy] = (ploidyval, purval, Aacc)
    return calldata

def readWGSEvidence(wgs):
    wgsdata = {}
    for line in open(wgs, "r"):
        if "Patient" in line:
            continue
        (patient, sample, tomcall, vafcat) = line.rstrip().split('\t')
        if patient not in wgsdata:
            wgsdata[patient] = {}
        wgsdata[patient][sample] = (tomcall, vafcat)
    return wgsdata

def readFinalCalls(callfile):
    finalcalls = {}
    for line in open(callfile, "r"):
        if "Patient" in line:
            continue
        (fullcall, patient, sample, NYGC) = line.rstrip().split('\t')
        if patient not in finalcalls:
            finalcalls[patient] = {}
        if "Nondiploid" in fullcall or "etraploid" in fullcall:
            simplecall = "Tetraploid"
        elif "iploid" in fullcall or "Unclear" in fullcall:
            simplecall = "Diploid"
        else:
            print("Unknown call", fullcall)
            assert(False)
        
        finalcalls[patient][sample] = (NYGC, simplecall)
    return finalcalls

def readGoodnessData():
    goodnesses = {}
    files = []
    for (__, __, f) in walk(goodness_dir + "500/diploid/"):
        files += f
    for f in files:
        if "goodness" not in f:
            continue
        patient = f.split("_")[0]
        for ploidy in ("diploid", "tetraploid"):
            if not isfile(goodness_dir + "500/" + ploidy + "/" + f):
                continue
            gfile = open(goodness_dir + "500/" + ploidy + "/" + f)
            if patient=="772":
                gfile = open(goodness_dir + "550/" + ploidy + "/" + f)
            for line in gfile:
                if "x" in line:
                    continue
                (pid, goodness) = line.rstrip().split()
                sample = pid.split('_')[1].split('"')[0]
                goodness = float(goodness)
                if patient not in goodnesses:
                    goodnesses[patient] = {}
                    goodnesses[patient]["diploid"] = {}
                    goodnesses[patient]["tetraploid"] = {}
                goodnesses[patient][ploidy][sample] = goodness
    return goodnesses


def getGoodnessDiff(patient, sample, goodnesses):
    diploid = 0
    tetraploid = 0
    if patient in goodnesses:
        if sample in goodnesses[patient]["diploid"]:
            diploid = goodnesses[patient]["diploid"][sample]
        if sample in goodnesses[patient]["tetraploid"]:
            tetraploid = goodnesses[patient]["tetraploid"][sample]
    if diploid==0 or tetraploid==0:
        return "Only one"
    if diploid - tetraploid > 0:
        return "Diploid"
    return "Tetraploid"

def getWGSEvidence(patient, sample, wgsdata):
    if patient in wgsdata and sample in wgsdata[patient]:
        return wgsdata[patient][sample]
    return ("Unknown", "Unknown")

def getDtMatches(patient, sample, flowdata, calldata):
    dmatch = "Unknown"
    tmatch = "Unknown"
    if patient in flowdata and patient in calldata:
        dmatch = "None"
        tmatch = "None"
        flowlist = flowdata[patient][2]
        assert(sample in calldata[patient])
        if "diploid" in calldata[patient][sample]:
            dploidy = calldata[patient][sample]["diploid"][0]
            dmatch = "False"
            #This is where we'd add in something about 'the 2N peak was pretty wide' if we got that information.
            if abs(dploidy-2) < 0.3:
                dmatch = "Two"
            else:
                for flow in flowlist:
                    if abs(flow-dploidy) < 0.2:
                        dmatch = "True"
                        break
        tploidy = "None"
        if "eight" in calldata[patient][sample]:
            tploidy = calldata[patient][sample]["eight"][0]
        elif "tetraploid" in calldata[patient][sample]:
            tploidy = calldata[patient][sample]["tetraploid"][0]
        if tploidy != "None":
            tmatch = "False"
            for flow in flowlist:
                if abs(flow-tploidy) < 0.2:
                    tmatch = "True"
                    break
    return (dmatch, tmatch)

def getBetterAccuracy(patient, sample, calldata):
    better_accuracy = "Unknown"
    if patient in calldata and sample in calldata[patient]:
        if "diploid" in calldata[patient][sample]:
            dacc = calldata[patient][sample]["diploid"][2]
            tacc = 0
            whichtet = "Tetraploid"
            if "tetraploid" in calldata[patient][sample]:
                tacc = calldata[patient][sample]["tetraploid"][2]
            if "eight" in calldata[patient]:
                whichtet = "Eight"
                tacc = calldata[patient][sample]["eight"][2]
            if tacc==0:
                better_accuracy = "Diploid only"
            else:
                if abs(dacc-tacc) < .03:
                    better_accuracy = "Neither"
                elif dacc > tacc:
                    better_accuracy = "Diploid"
                else:
                    better_accuracy = whichtet
        elif "eight" in calldata[patient][sample]:
            better_accuracy = "Eight only"
        elif "tetraploid" in calldata[patient][sample]:
            better_accuracy = "Tetraploid only"
    if patient=="772" and sample=="24571":
        #We had to go back and get a whole different gamma, since there wasn't a tet version, but we wanted one.
        return "Diploid only"
    return better_accuracy

def getAccuracyDiff(patient, sample, calldata):
    diff = 0.0
    if patient in calldata and sample in calldata[patient]:
        if "diploid" in calldata[patient][sample]:
            dacc = calldata[patient][sample]["diploid"][2]
            tacc = 0
            if "tetraploid" in calldata[patient][sample]:
                tacc = calldata[patient][sample]["tetraploid"][2]
            if "eight" in calldata[patient]:
                tacc = calldata[patient][sample]["eight"][2]
            diff = dacc-tacc
    return diff

def getFinalCall(patient, sample, finaldata):
    if patient in finaldata and sample in finaldata[patient]:
        return finaldata[patient][sample][1]
    return "Unknown"

def getFlowRatio(patient, flowdata):
    if patient in flowdata:
        return flowdata[patient][0]
    return "0::0"

def NYGCCloserTo(patient, sample, finaldata, calldata):
    NYGCcall = "NA"
    if patient in finaldata and sample in finaldata[patient]:
        NYGCcall = finaldata[patient][sample][0]
    if NYGCcall=="NA":
        return "NA"
    NYGCcall = float(NYGCcall)
    dploidy = "NA"
    tploidy = "NA"
    if patient in calldata and sample in calldata[patient]:
        if "diploid" in calldata[patient][sample]:
            dploidy = calldata[patient][sample]["diploid"][0]
        if "tetraploid" in calldata[patient][sample]:
            tploidy = calldata[patient][sample]["tetraploid"][0]
    if dploidy=="NA" and NYGCcall > 2.5:
        return "Tetraploid"
    if tploidy=="NA" and NYGCcall < 2.8:
        return "Diploid"
    if dploidy=="NA" or tploidy=="NA":
        assert(False)
    if abs(tploidy-NYGCcall) < abs(dploidy-NYGCcall):
        return "Tetraploid"
    else:
        return "Diploid"

def writeHeaders(outdata):
    outdata.write("Patient")
    outdata.write("\tSample")
    outdata.write("\tTom's Partek Call")
    outdata.write("\t2N VAF histogram category")
    outdata.write("\tFlow ratio")
    outdata.write("\tClose diploid flow?")
    outdata.write("\tClose tetraploid flow?")
    outdata.write("\tBetter accuracy?")
    outdata.write("\tAccuracy difference")
    outdata.write("\tNYGC closer to:")
    outdata.write("\tGoodness diff:")
    outdata.write("\tFinal Call")
    outdata.write("\n")
    

flowdata = readFlowSummary(flow_summary)
calldata = readCallSummary(call_summary)
wgsdata = readWGSEvidence(initcall)
finaldata = readFinalCalls(finalcall_file)
goodness_data = readGoodnessData()
outdata = open(outfile, "w")
writeHeaders(outdata)
for patient in calldata.keys():
    for sample in calldata[patient].keys():
        if not include_challenge and (patient not in wgsdata or sample not in wgsdata[patient]):
            continue
        (tomcall, vafcat) = getWGSEvidence(patient, sample, wgsdata)
        (dmatch, tmatch) = getDtMatches(patient, sample, flowdata, calldata)
        better_accuracy = getBetterAccuracy(patient, sample, calldata)
        flowratio = getFlowRatio(patient, flowdata)
        accuracy_diff = getAccuracyDiff(patient, sample, calldata)
        nygccloser = NYGCCloserTo(patient, sample, finaldata, calldata)
        goodness = getGoodnessDiff(patient, sample, goodness_data)
        final_call = getFinalCall(patient, sample, finaldata)
        writeNewLine(outdata, patient, (sample, tomcall, vafcat, flowratio, dmatch, tmatch, better_accuracy, accuracy_diff, nygccloser, goodness, final_call))
    
outdata.close()

