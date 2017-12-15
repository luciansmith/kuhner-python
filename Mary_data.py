# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 16:18:22 2016

@author: lpsmith
"""

import lucianSNPLibrary as lsl

prune = open("prunebreaks.txt", "r")
data = []
for line in prune:
    point = int(line)
    data.append(point)

print max(data)
binwidth = 1
label = "position"

lsl.createPrintAndSaveHistogram(data, "mary_out.txt", binwidth, xdata=label)
