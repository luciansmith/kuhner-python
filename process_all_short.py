# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 13:56:19 2016

@author: lpsmith
"""


def processThreeLines(threelines, output):
    line1 = threelines[0].rstrip().split("\t")
    line2 = threelines[1].rstrip().split("\t")
    line3 = threelines[2].rstrip().split("\t")
    if (line1[6] == "Double_d" and line2[6] == "Balanced_gain" and line3[6] == "Balanced_gain"):
        print line1, line2, line3
    key= line1[6] + "\t" + line2[6] + "\t" + line3[6]
    id = line1[1] + "_" + line1[2]
    if output.get(key) == None:
        output[key] = [id]
    else:
        output[key].append(id)

shortfile = open("short_segments/all_short.txt", "r")
threelines = []
output = {}
for line in shortfile:
    if (line.find("patient") != -1):
        continue
    if len(threelines) == 3:
        processThreeLines(threelines, output)
        threelines = []
    else:
        threelines.append(line)

for info in output.items():
    print info[0], "\t", len(info[1])