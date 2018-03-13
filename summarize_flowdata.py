# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 16:35:09 2018

@author: Lucian
"""
flowout = "flow_summary.txt"

flowfile = open("flow_raw_data.txt", "r")

flowdata = {}

for line in flowfile:
    if "Patient" in line:
        continue
    lvec = line.split()
    (patient, age, level, id, fourN) = lvec[0:5]
    if patient not in flowdata:
        flowdata[patient] =  [0, 0, set()]
    flowdata[patient][0] += 1
    aneuploidFlow = False
    fourN = float(fourN)
    if fourN > 6:
        aneuploidFlow = True
        flowdata[patient][2].add(4.0)
    if len(lvec) > 5:
        aneuploidFlow = True
        flowdata[patient][2].add(float(lvec[5]))
    if len(lvec) > 7:
        aneuploidFlow = True
        flowdata[patient][2].add(float(lvec[7]))
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