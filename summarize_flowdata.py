# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 16:35:09 2018

@author: Lucian
"""
flowout = "flow_summary.txt"

#flowfile = open("flow_raw_data.txt", "r")
flowfile = open("ChallengeFlowLucian_09202018_deidentified.txt", "r")

flowdata = {}

for line in flowfile:
    if "RandomID" in line:
        continue
    lvec = line.split('\t')
    (patient, age, flownum, sampnum, fourN) = lvec[0:5]
    if patient not in flowdata:
        flowdata[patient] =  [0, 0, set()]
    flowdata[patient][0] += 1
    aneuploidFlow = False
    if (fourN==""):
        fourN = 0
    fourN = float(fourN)
    if fourN > 6:
        aneuploidFlow = True
        flowdata[patient][2].add(4.0)
    (a1p, perca1, a2p, perca2) = lvec[5:9]
    if a1p != "":
        aneuploidFlow = True
        flowdata[patient][2].add(float(a1p))
    if a2p != "":
        aneuploidFlow = True
        flowdata[patient][2].add(float(a2p))
    if aneuploidFlow:
        flowdata[patient][1] += 1

flowfile = open(flowout, "w")
flowfile.write("Patient")
flowfile.write("\tnAneuploid")
flowfile.write("\tnTotal")
flowfile.write("\tSampled Aneuploids")
flowfile.write("\n")
for patient in flowdata:
    flowfile.write(patient)
    flowfile.write("\t" + str(flowdata[patient][1]) + ":" + str(flowdata[patient][0]))
    for val in flowdata[patient][2]:
        flowfile.write("\t" + str(val))
    flowfile.write("\n")
flowfile.close()