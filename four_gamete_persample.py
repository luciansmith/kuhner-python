#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:02:52 2016

@author: lpsmith
"""

#Create parallel BAF segmentation files to match the CN segmentation files

from __future__ import division
from os import walk
from os.path import isfile

import lucianSNPLibrary as lsl
import random
import numpy

#Use if you want to analyze the 'rejoined' data:
#expands_directory = "expands_rejoined_input/results_rejoined/"
#outDirectory = "expands_pairs_rejoined_shuffled/"

#Use if you want to analyze the original Xiaohong-segmented data:
expands_directory = "expands_input_shuffled_results/"
outDirectory = "expands_inshuffled_summary/"

mindiff = 10000
#deq_out = open("diseqs_test.txt", "w")
persample = open("four_gamete_inshuffled_persample.txt", "w")
perpatient = open("four_gamete_inshuffled_perpatient.txt", "w")

def combinePopsInList(combinedpops, point, sample, pops):
    if not(point in combinedpops):
        combinedpops[point] = {}
    combinedpops[point][sample] = pops

def findFourGametes(segments):
    foundfour = 0
    notfour = 0
    for n in range(0, len(segments)-1):
        gametes1 = segments[n][3]
        for s in range(n+1, len(segments)):
            gametes2 = segments[s][3]
            if (len(gametes1) != len(gametes2)):
                print "Error!"
            g4test = set()
            g4test.add("zerozero") #Because we know that the ancestral allele was WT.
            for g in range(0,len(gametes1)):
                g1 = gametes1[g]
                g2 = gametes2[g]
                if (g1==1):
                    if (g2==1):
                        g4test.add("oneone")
                    else:
                        g4test.add("onezero")
                else:
                    if (g2==1):
                        g4test.add("zeroone")
                    else:
                        g4test.add("zerozero")
            if len(g4test) > 3:
                foundfour += 1
            else:
                notfour += 1
    return (foundfour, notfour)


def findFourGametesPopsOnly(allpops):
    foundfour = 0
    notfour = 0
    for n in range(0, len(allpops)-1):
        gametes1 = allpops[n]
        for s in range(n+1, len(allpops)):
            gametes2 = allpops[s]
            if (len(gametes1) != len(gametes2)):
                print "Error!"
            g4test = set()
            g4test.add("zerozero") #Because we know that the ancestral allele was WT.
            for g in range(0,len(gametes1)):
                g1 = gametes1[g]
                g2 = gametes2[g]
                if (g1==1):
                    if (g2==1):
                        g4test.add("oneone")
                    else:
                        g4test.add("onezero")
                else:
                    if (g2==1):
                        g4test.add("zeroone")
                    else:
                        g4test.add("zerozero")
            if len(g4test) > 3:
                foundfour += 1
            else:
                notfour += 1
    return (foundfour, notfour)


flist = []
filesets = {}
#SNPfiles.append(["1034", "20008"])
for (_, _, f) in walk(expands_directory):
    flist += f
    
for f in flist:
    if f.find(".sps") != -1 and f.find(".cbs") == -1 and f.find(".spstats") == -1:
        treename = f.replace(".sps", ".tree")
        if not(isfile(expands_directory + treename)):
            continue
        (patient, sample, tag) = f.split("_")
        if not(patient in filesets):
            filesets[patient] = list()
        filesets[patient].append(f)

persample.write("Patient\tSample\tnumpops\tnFourGametes\tnNotFour\tFoundRandomFreq\tavgFourGametes\tstDevFourGametes\tDiseq\tavg_diseq\tstd_diseq\n")
patientdata = {}
for patient in filesets:
    fileset = filesets[patient]
    if len(fileset) == 1:
        continue
    numsamples = 0
    for outfile in fileset:
        (patient, sample, __) = outfile.split("_")
        efile = open(expands_directory + outfile, "r")
        numpops = 0
        segments = []
        for line in efile:
            if line.find("expands version") != -1:
                continue
            if line.find("chr") != -1:
                numpops = len(line.split()) - 16
                continue
            rowvec = line.rstrip().split()
            scenario = rowvec[14]
            if (scenario == "4"):
                #Don't count the ones that separate the CN and the BAF segments into different populations.
                continue
            chr = int(rowvec[0])
            start = int(float(rowvec[1]))
            end = int(float(rowvec[2]))
            pops = rowvec[15:15+numpops]
            for p in range(0,len(pops)):
                pops[p] = int(pops[p])
            segment = [chr, start, end, pops]
            segments.append(segment)
        if (numpops <3):
            perpatient.write(patient + "\t" + sample + "\t" + str(numpops) + "\t---\n")
            continue
        (foundfour, notfour) = findFourGametes(segments)
        numlteqFF = 0
        numtot = 10000
        totfoundfour = []
        samplepops = []
        for segment in segments:
            samplepops.append(segment[3])
        samplepops = numpy.array(samplepops)
        for n in range(0,numtot):
            lsl.shufflePops(segments)
            (rfoundfour, rnotfour) = findFourGametes(segments)
            if (rfoundfour <= foundfour):
                numlteqFF += 1
            totfoundfour.append(rfoundfour)
        print "Calculating info for patient", patient, "sample", sample
        probones = []
        numsegs = len(samplepops)
        numpops = len(samplepops[0])
        for n in range(0, numsegs):
            probones.append(numpy.sum(samplepops[n])/len(samplepops[n]))
        persample.write(patient + "\t" + sample + "\t" + str(numpops) + "\t" + str(foundfour) + "\t" + str(notfour) + "\t" + str(numlteqFF/numtot) + "\t" + str(numpy.average(totfoundfour)) + "\t" + str(numpy.std(totfoundfour)) + "\n")
#        avgdiseq = lsl.CalculateAverageDiseq(samplepops, probones, numsegs, numpops)
#        shuffledaverages = []
#        shuffled = samplepops
#        for n in range(0,1000):
#            shuffled = lsl.ShufflePops(shuffled)
#            avg = lsl.CalculateAverageDiseq(shuffled, probones, numsegs, numpops)
#            shuffledaverages.append(avg)
#        print "Writing info for patient", patient, "sample", sample
#        persample.write(patient + "\t" + sample + "\t" + str(numpops) + "\t" + str(foundfour) + "\t" + str(notfour) + "\t" + str(numlteqFF/numtot) + "\t" + str(numpy.average(totfoundfour)) + "\t" + str(numpy.std(totfoundfour)) + "\t" + str(avgdiseq) + "\t" + str(numpy.average(shuffledaverages)) + "\t" + str(numpy.std(shuffledaverages)) + "\n")
        if not(patient in patientdata):
            patientdata[patient] = []
        ss = (sample, segments)
        patientdata[patient].append(ss)
        numsamples += 1
    armlist = []
    if (patient in patientdata):
        armlist = lsl.OneSegmentPerArm(patientdata[patient], numsamples, mindiff)
    if (len(armlist) < 2):
        perpatient.write(str(patient) + "\t" + str(numsamples) + "\t" + str(len(armlist)) + "\t---\t---\t---\n")
        continue
    allpops = []
    for seglist in armlist:
        onearmpops = []
        for oneseg in seglist[1]:
            onearmpops += oneseg[1][3]
        allpops.append(onearmpops)
    allpops = numpy.array(allpops)
    (foundfour, notfour) = findFourGametesPopsOnly(allpops)
    numlteqFF = 0
    numtot = 10000
    totfoundfour = []
    samplepops = []
    for segment in allpops:
        samplepops.append(segment[3])
    samplepops = numpy.array(samplepops)
    for n in range(0,numtot):
        randomizePopulationsPopsOnly(allpops)
        (rfoundfour, rnotfour) = findFourGametesPopsOnly(allpops)
        if (rfoundfour <= foundfour):
            numlteqFF += 1
        totfoundfour.append(rfoundfour)
    print "Calculating info for patient", patient, "sample", sample
    probones = []
    numsegs = len(samplepops)
    numpops = len(samplepops[0])
    for n in range(0, numsegs):
        probones.append(numpy.sum(samplepops[n])/len(samplepops[n]))
    perpatient.write(patient + "\t" + str(numpops) + "\t" + str(foundfour) + "\t" + str(notfour) + "\t" + str(numlteqFF/numtot) + "\t" + str(numpy.average(totfoundfour)) + "\t" + str(numpy.std(totfoundfour)) + "\n")
    #print allpops
#    probones = []
#    numsegs = len(allpops)
#    numpops = len(allpops[0])
#    for n in range(0, numsegs):
#        probones.append(numpy.sum(allpops[n])/len(allpops[n]))
#    avgdiseq = lsl.CalculateAverageDiseq(allpops, probones, numsegs, numpops)
#    shuffledaverages = []
#    shuffled = allpops
#    for n in range(0,1000):
#        shuffled = lsl.ShufflePops(shuffled)
#        avg = lsl.CalculateAverageDiseq(shuffled, probones, numsegs, numpops)
#        shuffledaverages.append(avg)
#        deq_out.write(str(avg))
#        deq_out.write("\t")
#    deq_out.write("\n")
#    perpatient.write(str(patient) + "\t" + str(numsamples) + "\t" + str(len(allpops)) + "\t" + str(avgdiseq) + "\t" + str(numpy.average(shuffledaverages)) + "\t" + str(numpy.std(shuffledaverages)) + "\n")

#deq_out.close()
perpatient.close()
persample.close()

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        