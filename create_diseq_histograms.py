# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:33:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk

import lucianSNPLibrary as lsl
import numpy

# read the filtered data that compares Xiaohong's segmentation data with raw SNP data

#filenames = ["1049_20780_avglog2rs.txt", "1049_20782_avglog2rs.txt"]
filename = "diseqs.txt"
all_data = []



file = open(filename, "r")
id = 0
all_data = []
for line in file:
    id += 1
    all_data = numpy.array(map(float, line.rstrip().split()))
    #binwidth = (max(all_data) - min(all_data))/100
    #binwidth = pow(10,int(numpy.floor(numpy.log10(abs(binwidth)))))
    binwidth = 0.001
    lsl.createPrintAndSaveHistogram(all_data, "diseq_" + str(id) + ".txt", binwidth, xdata="diseq")
