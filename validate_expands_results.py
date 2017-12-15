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
from os.path import isdir
import sys

import lucianSNPLibrary as lsl
import numpy

#Use if you want to analyze the 'rejoined' data:
#expands_directory = "expands_rejoined_input/results_rejoined/"

#Use if you want to analyze the original Xiaohong-segmented data:
#expands_directory = "expands_shuffled_results/"

#Use for full expands run, non-rejoined
expands_core = "shuffled0/results/"

#Use for full expands run, rejoined
#expands_core = "expands_full_rejoined_results"

if len(sys.argv) == 2:
    expands_core = sys.argv[1]
    
elif len(sys.argv) > 2:
    print sys.argv[0], "can only take one argument (or none, if you want to use the default 'expands_full_results')."
    exit(1)

if not(isdir(expands_core)):
    print "'" + expands_core + "' is not a valid directory."
    exit(1)

expands_directory = expands_core
if expands_directory[len(expands_directory)-1] != "/":
    expands_directory = expands_directory + "/"

slash = expands_core.find("/")
if (slash != -1):
    expands_core = expands_core[0:slash]


mindiff = 10000
numtot = 500 #Number of times to shuffle the populations
#deq_out = open("diseqs_test.txt", "w")
persample = open("validate_" + expands_core + "_persample.txt", "w")
perpatient = open("validate_" + expands_core + "_perpatient.txt", "w")



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
#        if (patient != "1072"):
#            continue
        if not(patient in filesets):
            filesets[patient] = list()
        filesets[patient].append(f)

patientdata = {}
persample.write("Patient\tSample\tnumpops\tnumsegs\tnFourGametes\tnNotFour\tFoundRandomFreq\tavgFourGametes\tstdFourGametes\n")
perpatient.write("Patient\tnumsamples\tnumsegs\tnumpops\tnFourGametes\tnNotFour\tFoundRandomFreq\tavgFourGametes\tstdFourGametes\n")
print "Patient\tSample(s)\tnumpops\tnumsegs\tnFourGametes\tnNotFour\tFoundRandomFreq\tavgFourGametes\tstdFourGametes"
for patient in filesets:
    fileset = filesets[patient]
    if len(fileset) == 1:
        continue
    numsamples = 0
    for outfile in fileset:
        (patient, sample, __) = outfile.split("_")
        efile = open(expands_directory + outfile, "r")
        #print outfile
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
        justpops = lsl.GetJustPops(segments)
        (foundfour, notfour) = lsl.FindFourGametes(justpops)
        numlteqFF = 0
        totfoundfour = []
        samplepops = []
        for segment in segments:
            samplepops.append(segment[3])
        samplepops = numpy.array(samplepops)
        shuffled = justpops
        (goodpops, goodsegs) = lsl.GetVariedRowsAndColumns(shuffled)
        if (len(goodpops)==0 or len(goodsegs)==0):
            numlteqFF = numtot
            totfoundfour = [foundfour]
        else:
            for n in range(0,numtot):
                #if n/10 == numpy.floor(n/10):
                    #    print "shuffling", n
                    shuffled = lsl.ShufflePops(shuffled, goodpops, goodsegs)
                    (rfoundfour, rnotfour) = lsl.FindFourGametes(shuffled)
                    if (rfoundfour <= foundfour):
                        numlteqFF += 1
                    totfoundfour.append(rfoundfour)
        #print "Calculating info for patient", patient, "sample", sample
        probones = []
        numsegs = len(samplepops)
        numpops = len(samplepops[0])
        for n in range(0, numsegs):
            probones.append(numpy.sum(samplepops[n])/len(samplepops[n]))
        persample.write(patient + "\t" + sample + "\t" + str(numpops) + "\t" + str(numsegs) + "\t" + str(foundfour) + "\t" + str(notfour) + "\t" + str(numlteqFF/numtot) + "\t" + str(numpy.average(totfoundfour)) + "\t" + str(numpy.std(totfoundfour)) + "\n")
        print patient + "\t" + sample + "\t" + str(numpops) + "\t" + str(numsegs) + "\t" + str(foundfour) + "\t" + str(notfour) + "\t" + str(numlteqFF/numtot) + "\t" + str(numpy.average(totfoundfour)) + "\t" + str(numpy.std(totfoundfour))
        
        if not(patient in patientdata):
            patientdata[patient] = []
        ss = (sample, segments)
        patientdata[patient].append(ss)
        numsamples += 1
    armlist = lsl.OneSegmentPerArm(patientdata[patient], numsamples, mindiff)
    if (len(armlist) < 2):
        perpatient.write(str(patient) + "\t" + str(numsamples) + "\t" + str(len(armlist)) + "\t---\t---\t---\n")
        print str(patient) + "\t" + str(numsamples) + "\t" + str(len(armlist)) + "\t---\t---\t---"
        continue
    allpops = []
    for seglist in armlist:
        onearmpops = []
        for oneseg in seglist[1]:
            onearmpops += oneseg[1][3]
        allpops.append(onearmpops)
    allpops = numpy.array(allpops)
    (foundfour, notfour) = lsl.FindFourGametes(allpops)
    #print allpops
    probones = []
    totfoundfour = []
    numsegs = len(allpops)
    numpops = len(allpops[0])
    for n in range(0, numsegs):
        probones.append(numpy.sum(allpops[n])/len(allpops[n]))
    avgdiseq = lsl.CalculateAverageDiseq(allpops, probones, numsegs, numpops)
    shuffledaverages = []
    shuffled = allpops
    (goodpops, goodsegs) = lsl.GetVariedRowsAndColumns(shuffled)
    numlteqFF = 0
    for n in range(0,numtot):
        shuffled = lsl.ShufflePops(shuffled, goodpops, goodsegs)
        (rfoundfour, rnotfour) = lsl.FindFourGametes(shuffled)
        if (rfoundfour <= foundfour):
            numlteqFF += 1
        totfoundfour.append(rfoundfour)
        
    
#    print str(patient) + "\t" + str(numsamples) + "\t" + str(len(allpops)) + "\t" + str(avgdiseq) + "\t" + str(numpy.average(shuffledaverages)) + "\t" + str(numpy.std(shuffledaverages))
    perpatient.write(patient + "\t" + str(len(patientdata[patient])) + "\t" + str(numpops) + "\t" + str(len(armlist)) + "\t" + str(foundfour) + "\t" + str(notfour) + "\t" + str(numlteqFF/numtot) + "\t" + str(numpy.average(totfoundfour)) + "\t" + str(numpy.std(totfoundfour)) + "\n")
    print patient + "\t" + str(len(patientdata[patient])) + "\t" + str(numpops) + "\t" + str(len(armlist)) + "\t" + str(foundfour) + "\t" + str(notfour) + "\t" + str(numlteqFF/numtot) + "\t" + str(numpy.average(totfoundfour)) + "\t" + str(numpy.std(totfoundfour))
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
persample.close()
perpatient.close()
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        