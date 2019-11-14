#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 12:45:13 2019

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages/pymix")

import imp
drp = imp.load_source("drp","/home/mkkuhner/Papers/phylo/dropout.py")

import numpy
import operator
import mixture

import lucianSNPLibrary as lsl

onlysomepatients = True
somepatients = ["521"]

VAFdir = "VAFclusters/"
outdir = "mary_histpeaks/"
#outdir = "VAFclusters_histograms/"
#outdir_low = "VAFclusters_histograms_low_oneplus_521/"

###########################################################################
class peak:
  def  __init__(self, lucianmax, histmax, histsites):
    self.lucianmax = lucianmax
    self.histmax = histmax
    self.histsites = histsites

# label is a string; get tips out of it
def split_label(label):
  # remove parentheses, single quotes
  newlabel = ""
  for char in label:
    if char not in ["(",")","'"," "]:
      newlabel = newlabel + char
  # split on commas
  label = newlabel.split(",")
  items = []
  for item in label:
    if item == "":  continue
    items.append(item)
  return tuple(items)

def dropchance(vaf,totreads):
  return (1.0 - vaf)**totreads

def sortLabels(labels):
    newlist = []
    sublist = []
    for label in labels:
        if len(label.split(", '"))==1:
            sublist.append(label)
    sublist.sort()
    newlist.extend(sublist)
    for n in range(9,1,-1):
#    for n in range(2,9):
        sublist = []
        for label in labels:
            if len(label.split(", '"))==n:
                sublist.append(label)
        sublist.sort()
        newlist.extend(sublist)
    return newlist

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

def makeFilename(label):
    elements = eval(label)
    ret = ""
    for element in elements:
        ret += element + "-"
    return ret[:-1]

def getHistMaxes(hist,verbose):
    ret = []
    keylist = list(hist.keys())
    keylist.sort()
    localmax = 0
    maxkey = keylist[0]
    localmin = 0
    direction = "up"
    distance = 0
    for key in keylist:
        val = hist[key]
        if direction=="up":
            if val > localmax:
                localmax = val
                maxkey = key
                distance = 0
            else:
                distance += 1
            if distance >= 40 and localmax-val > 0.25:
                if verbose:  print("Switched directions: going down at", key, distance, localmax, val)
                direction = "down"
                if verbose:  print("Found value at",maxkey,localmax)
                ret.append(maxkey)
                localmin = val
        elif direction=="down":
            if val < localmin:
                localmin = val
                distance = 0
            else:
                distance += 1
            if distance >= 40 and val - localmin > 0.25:
                if verbose:  print("Switched directions: going up at", key, distance, localmax, val)
                direction="up"
                localmax = val
                maxkey = key
    return ret

def getPymixFail(libmaxes, histpeaks):
#                    summary.write("\t"+str(libmaxes[i][0]))
#                    summary.write("\t"+str(libmaxes[i][1]))
    libpeaks = []
    for (peak, height) in libmaxes:
        for oldpeak in libpeaks:
            if abs(oldpeak-peak) < 0.01:
                print("Failure: two libpeaks are too close together.", str(oldpeak), str(peak))
                return True
        libpeaks.append(peak)
    for peak in libpeaks:
        foundclose = False
        for histpeak in histpeaks:
            if abs(histpeak-peak) < 0.01:
                foundclose = True
        if not foundclose:
            print("Failure: no histogram peak found close enough to", str(peak), "in", str(histpeaks))
            return True
    return False


##############################################################################

import os
if not os.path.exists(outdir):
    os.mkdir(outdir)
#if not path.isdir(outdir_low):
#    mkdir(outdir_low)

VAFfiles = []
for __, _, files in walk(VAFdir):
    VAFfiles += files

(patientSampleMap, samplePatientMap) = lsl.getPatientSampleMap(dipvtet_file="calling_evidence_odds.tsv")
deletions, CNVs = lsl.loadDeletionsAndCNVs(samplePatientMap)

summary = open("summary_smoothed_and_fit.tsv", "w")
summary.write("Patient")
summary.write("\tSample")
summary.write("\tnPoints")
summary.write("\tCall")
summary.write("\tGroup")
#summary.write("\tMean x2")
#summary.write("\tStdev x2")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\tHistMax x2")
#summary.write("\tHistMax height")
summary.write("\tFitNormal x2")
summary.write("\tFitNormal weight")
summary.write("\n")

# new data structure to make Mary's tables
# dictionary by patient, sample, partition, call, contains peaks
mdata = {}

# postdict is a dictionary [patient][chr][pos] of posrecords (one per
# mutant position) with each posrecord eventually to contain mutrecords--
# but not yet as they rely on information not in hand, so we make them
# empty.
posdict = {}

for file in VAFfiles:
    partitions = []
    if "_VAFs" not in file:
        continue
    (patient, sample) = file.split("_")[0:2]
    if onlysomepatients and patient not in somepatients:
        continue
    if patient not in mdata:
      mdata[patient] = {}
      posdict[patient] = {}
    if sample not in mdata[patient]:
      mdata[patient][sample] = {}
    data = {}
    for line in open(VAFdir + file, "r"):
        lvec = line.rstrip().split("\t")
        if "Patient" in line:
            labels = lvec[5:]
            for group in labels:
                data[group] = {}
            continue
        if lvec[2]=="23" or lvec[2]=="24":
            continue
        chrom = lvec[2]
        pos = int(lvec[3])

        if chrom not in posdict[patient]:
          posdict[patient][chrom] = {}

        call = getCNVCall(patient, sample, lvec[2], int(lvec[3]), CNVs)
        if call==(-1, -1):
            continue
        for n in range(5, len(lvec)):
            label = labels[n-5]
            if lvec[n] != "":
               # debug
               if pos not in posdict[patient][chrom]:
                 posdict[patient][chrom][pos] = drp.posrecord(split_label(label))
               if call not in data[label]:
                    data[label][call] = []
               data[label][call].append(float(lvec[n]))
    for label in data:
        if label not in mdata[patient][sample]:
            mdata[patient][sample][label] = {}
        for call in data[label]:
            if len(data[label][call]) == 0:
                continue
            if len(data[label][call]) < 100:
                #Skip groups with fewer than 100 VAFs.
               continue
            if call not in mdata[patient][sample][label]:
               mdata[patient][sample][label][call] = []
            filename = patient + "_" + sample + "_" + makeFilename(label) + "_" + str(call[0]) + "_" + str(call[1]) + "_hist.png"
            # DEBUG
            group = split_label(label)
            if len(group) == 4 and sample == "23579":
              verbose = True
            else:  verbose = False

            if not verbose:
              hist = lsl.createPrintAndSaveHistogram(data[label][call], filename, 0.001, xdata="VAF", savefig=False, show=False)
            else: 
              hist = lsl.createPrintAndSaveHistogram(data[label][call], filename, 0.001, xdata="VAF", savefig=False, show=True)
            mean = numpy.mean(data[label][call])
            stdev = numpy.std(data[label][call])
            histmaxes = getHistMaxes(hist,verbose)
            print("Hitmaxes",histmaxes)
            #print(patient, sample, label, call)
            ###THIS IS WHERE YOU FIND THE HISTOGRAM PEAKS###
            ##Data:  data[label][call]
            ##Peaks:  histmaxes
            ##Peak heights:  hist[histmaxes[n]]
            ##Stdev:  stdev

            emdata = mixture.DataSet()           
            emdata.fromList(data[label][call])
            numpeaks = len(histmaxes)
            gaussian_objects = []
            weights = []
            for i in xrange(numpeaks):
              n = mixture.NormalDistribution(histmaxes[i],stdev)
              gaussian_objects.append(n)
              weights.append(hist[histmaxes[i]])
            totweight = float(sum(weights))
            weights = [x/totweight for x in weights]
            mymix = mixture.MixtureModel(numpeaks,weights,gaussian_objects)
            mymix.EM(emdata,40,0.1)

            summary.write(patient)
            summary.write("\t" + sample)
            summary.write("\t" + str(len(data[label][call])))
            summary.write("\t" + str(call))
            summary.write("\t" + label)
            histmaxes.sort(reverse=True)
            libmaxes = []
            for i in range(mymix.G):
              libmax = 2*mymix.components[i].distList[0].mu
              lib_nsnv = mymix.pi[i]*len(data[label][call])
              libmaxes.append([libmax,lib_nsnv])
            libmaxes.sort(reverse=True)
            mydata = [patient,sample,label,call,len(data[label][call])]
            fail = getPymixFail(libmaxes, histmaxes)
            histHeightTotal = 0
            for hmax in histmaxes:
                histHeightTotal += hist[hmax]
            for i in range(mymix.G):
                histmax = histmaxes[i]
                summary.write("\t" + str(histmax*2))
                if (fail):
                    summary.write("\t[pymix failure]")
                    summary.write("\t"+str(len(data[label][call] * hist[histmax]/histHeightTotal)))
                else:
                    summary.write("\t"+str(libmaxes[i][0]))
                    summary.write("\t"+str(libmaxes[i][1]))
                thispeak = peak(histmax*2,libmaxes[i][0],libmaxes[i][1])
                mydata.append(thispeak)
            summary.write("\n")
            mdata[patient][sample][label][call] = mydata

# kanika data is a dictionary [pid][sid][chrom][pos][ref/alt] = count

# read stored kanika data
infile = open("kanika_readcounts.tsv","r")
filelines = infile.readlines()

hdrline = filelines[0].rstrip().split("\t")
pid_ind = hdrline.index("Patient")
sample_ind = hdrline.index("Sample")
chrom_ind = hdrline.index("chrom")
pos_ind = hdrline.index("pos")
alt_ind = hdrline.index("alt")
ref_ind = hdrline.index("ref")

kanika = {}
for line in filelines[1:]:
  line = line.rstrip().split("\t")
  pid = line[pid_ind]
  sample = line[sample_ind]
  chrom = line[chrom_ind]
  pos = int(line[pos_ind])
  alt = int(line[alt_ind])
  ref = int(line[ref_ind])

  if pid not in kanika:
    kanika[pid] = {}
  if sample not in kanika[pid]:
    kanika[pid][sample] = {}
  if chrom not in kanika[pid][sample]:
    kanika[pid][sample][chrom] = {}
  if pos not in kanika[pid][sample][chrom]:
    kanika[pid][sample][chrom][pos] = {}
  kanika[pid][sample][chrom][pos]["ref"] = ref
  kanika[pid][sample][chrom][pos]["alt"] = alt

# posdict is a dictionary [patient][chr][pos] of posrecords (one per
# mutant position) with each posrecord containing mutrecords for each
# sample.

for patient in kanika:
  if onlysomepatients and patient not in somepatients:  continue
  for sample in kanika[patient]:
    for chrom in kanika[patient][sample]:
      for pos in kanika[patient][sample][chrom]:
        pos = int(pos)
        k = kanika[patient][sample][chrom][pos]
        call = getCNVCall(patient, sample, chrom, pos, CNVs)
        if call == (-1,-1):  continue   # we can do nothing with call-unknown data
        try:
          p = posdict[patient][chrom][pos]
        except:
          #debu DEBUG
          print("failed to find")
          print(patient,chrom,pos)
          print("from kanika sample",sample)
          print("in posdict")
          exit()
        ref = k["ref"]
        alt = k["alt"]
        totreads = ref + alt
        vaf = alt/float(totreads)
        mut = drp.mutrecord(sample,call,vaf,totreads)
        p.addmut(sample,mut)

# sort by patients
patlist = mdata.keys()
patlist.sort()
wanted_call = (1,1)
partnames = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N"]
codes = ["a","b","c","d","e","f","g","h","i"]

for patient in patlist:
  # create codes for each sample in this patient
  samplelist = mdata[patient].keys()
  samplelist.sort()
  samplecodes = {}   # dictionary of sample ID to sample code
  for sample,code in zip(samplelist,codes):
    samplecodes[sample] = code
  
  # get ivafs, which is a dictionary[label][sample][call] containing
  # [list of vafs, list of counts]
  # made from mdata, which is a dictionary by patient, sample, label, 
  # call, and contains peak objects

  ivafs = {}
  for sample in mdata[patient]:
    for label in mdata[patient][sample]:
      if label not in ivafs:
        ivafs[label] = {}
      if sample not in ivafs[label]:
        ivafs[label][sample] = {}
      for call in mdata[patient][sample][label]:
        peaks = mdata[patient][sample][label][call][5:]
        allvafs = []
        allcounts = []
        for peak in peaks:
          allvafs.append(peak.lucianmax)
          allcounts.append(peak.histsites)
        ivafs[label][sample][call] = [allvafs,allcounts]

  # write the report
  # format:  header line with patient ID and partition codes and counts
  # one body line per sample, giving code, sid, VAF(s) and count(s),
  # with dropout numbers in square brackets

  outfilename = outdir+str(patient)+"_report.tsv"
  outfile = open(outfilename,"w")
  samplelist = mdata[patient].keys()
  samplelist.sort()
  # find all the non-deletion labels 
  all_labels = set()
  for sample in samplelist:
    for label in mdata[patient][sample]:
    # only partitions without a deletion entry
      if "-" in label:  continue
      for call in mdata[patient][sample][label]:
        all_labels.add(label)
  all_labels = list(all_labels)

  # sort the partitions into a sensible order; we'll try number of sites
  # in the partition in the first sample, summed over all calls, in
  # descending order
  partinfo = []

  for label in all_labels:
    group = split_label(label)
    totsites = 0
    code = [] 
    for entry in group:
      code.append(samplecodes[entry])
    code.sort()
    code = "".join(code)
    for call in mdata[patient][group[0]][label]:
      totsites += mdata[patient][group[0]][label][call][4]
    partinfo.append([label,group,code,totsites])
  partinfo = sorted(partinfo, key=lambda x: x[3], reverse=True)

  # write report header
  outline = "Patient" + patient + "\tSample"
  for entry in partinfo:
    outline += "\t" + entry[2] + "\tnMuts:" + entry[2]
  outline += "\n"
  outfile.write(outline)
  
  # find out what calls are relevant for each sample
  calls_by_sample = {}
  for sample in mdata[patient]:
    mycalls = set()
    for label in mdata[patient][sample]:
      for call in mdata[patient][sample][label]:
        mycalls.add(call)
    mycalls = list(mycalls)
    mycalls.sort()
    calls_by_sample[sample] = mycalls

  for sample in samplelist:
     for call in calls_by_sample[sample]:
       outline = samplecodes[sample] + "-" + str(call) + "\t" + sample 
       for label,group,code,totsites in partinfo:
         if sample not in ivafs[label]:
           outline += "\t--\t--"
           continue
         if call not in ivafs[label][sample]:
           outline += "\t--\t--"
           continue
         vafs, counts = ivafs[label][sample][call]
         outline += "\t" + str("%.2f" % round(vafs[0],2))
         for vaf in vafs[1:]:
           outline += ":" + str("%.2f" % round(vaf,2))
         outline += "\t" + str(int(counts[0]))
         for count in counts[1:]:
           outline += ":" + str(int(count))
       outline += "\n"
       outfile.write(outline)

  exit()

  # stopped here!

  # make map of partitions giving vafs for each
  vafmap = {}
  countmap = {}
  for label,group in zip(all_labels,all_groups):
    groupvafs = []
    groupcounts = []
    for sample in samplelist:
      if sample in group:
        peaks = mdata[patient][sample][label][wanted_call][5:]
        peakvafs = []
        peakcounts = []
        for peak in peaks:
          peakvafs.append(peak.lucianmax)
          peakcounts.append(peak.histsites)
      else:
        peakvafs = [0.0,]
        peakcounts = [0,]
      groupvafs.append(peakvafs)
      groupcounts.append(peakcounts)
    vafmap[group] = groupvafs
    countmap[group] = groupcounts
  
  for label,group,name in zip(all_labels,all_groups,partnames):
    outline = name + "\t" + label
    groupvafs = []
    for sample in samplelist:
      samplevafs = []
      mypeaks = []
      if sample in group:
        outline += "\t"
        peaks = mdata[patient][sample][label][wanted_call][5:]
        numpeaks = len(peaks)
        for i in xrange(numpeaks):
          outline += str(peaks[i].lucianmax)
          if i < numpeaks - 1:
            outline += ","
          mypeaks.append(peaks[i].lucianmax)
      else:
        outline += "\t---"
        mypeaks.append([0.0,])

    outline += "\n"
    outfile.write(outline)

  outfile.close()

