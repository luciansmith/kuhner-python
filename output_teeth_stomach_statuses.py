#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 12:10:15 2019

@author: lpsmith
"""

from __future__ import division
import csv

statusfile = "TS_statuses.tsv"

def getSampleStatuses():
    statuses = {}
    levels = {}
    t1t2_levels = {}
    with open("P01CA91955-WGS80-Full-Pilot-Samples.csv", "r") as csvfile:
        for lvec in csv.reader(csvfile):
            if "RandomID" in lvec:
                continue
            sample = lvec[6]
            statuses[sample] = {}
            statuses[sample]["patient"] = lvec[0]
            patient = lvec[0]
            statuses[sample]["prog"] = lvec[1]
            statuses[sample]["sex"] = lvec[2]
            statuses[sample]["age"] = lvec[5]
            time = lvec[8]
            if time=="LongFollowUp":
                time = "T3"
            if time=="Index":
                time = "T2"
            if sample == "23521":
                time = "T2"
            statuses[sample]["time"] = time
            level = lvec[7]
            gej = lvec[10]
            if gej != "NA":
                gej = int(gej)
            gej_dist = "NA"
            if level != "NA" and time in ["T1", "T2", "T3", "Pilot"]:
                level = int(level)
                gej_dist = gej-level
                if gej_dist<0:
                    gej_dist=0
                if patient not in levels:
                    levels[patient] = []
                    t1t2_levels[patient] = []
                levels[patient].append(gej_dist)
                if time in ["T1", "T2"]:
                    t1t2_levels[patient].append(gej_dist)
            statuses[sample]["level"] = level
            statuses[sample]["gej_dist"] = gej_dist
    for patient in levels:
        levels[patient].sort()
        t1t2_levels[patient].sort()
    for sample in statuses:
        patient = statuses[sample]["patient"]
        dist = statuses[sample]["gej_dist"]
        lcode = "S" #'stomach'
        t1t2_lcode = "S"
        assert(len(t1t2_levels[patient]) == 4)
        halfindex = round(len(levels[patient])/2)
        #print(str(halfindex))
        if dist=="NA":
            lcode = "NA"
            t1t2_lcode = "NA"
        else:
            if dist - levels[patient][0] > levels[patient][-1] - dist:
                lcode = "T" #'teeth'
            if dist - levels[patient][0] == levels[patient][-1] - dist:
                #print("Equal distances:", dist, levels[patient])
                if dist==levels[patient][halfindex]:
                    lcode = "T"
            if dist - t1t2_levels[patient][0] > t1t2_levels[patient][-1] - dist:
                t1t2_lcode = "T" #'teeth'
            if dist - t1t2_levels[patient][0] == t1t2_levels[patient][-1] - dist:
                #print("Equal distances:", dist, levels[patient])
                if dist==t1t2_levels[patient][2]:
                    t1t2_lcode = "T"
        statuses[sample]["levelcode_all"] = lcode
        statuses[sample]["levelcode_t1t2"] = t1t2_lcode
    return statuses


statuses = getSampleStatuses()

outfile = open(statusfile, "w")
outfile.write("Patient")
outfile.write("\tSample")
outfile.write("\tProg")
outfile.write("\tSex")
outfile.write("\tAge")
outfile.write("\tTime")
outfile.write("\tLevel")
outfile.write("\tGEJ_distance")
outfile.write("\tLevelCode_all")
outfile.write("\tLevelCode_T1T2")
outfile.write("\n")

for sample in statuses:
    outfile.write(statuses[sample]["patient"])
    outfile.write("\t" + sample)
    outfile.write("\t" + statuses[sample]["prog"])
    outfile.write("\t" + statuses[sample]["sex"])
    outfile.write("\t" + statuses[sample]["age"])
    outfile.write("\t" + statuses[sample]["time"])
    outfile.write("\t" + str(statuses[sample]["level"]))
    outfile.write("\t" + str(statuses[sample]["gej_dist"]))
    outfile.write("\t" + str(statuses[sample]["levelcode_all"]))
    outfile.write("\t" + str(statuses[sample]["levelcode_t1t2"]))
    outfile.write("\n")

outfile.close()