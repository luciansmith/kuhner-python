# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:24:20 2016

@author: lpsmith
"""

import numpy

#calculates the median BAF from array data.

midcutoffarray = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
allarraydata = []
midarraydata = [[]]
for c in range(len(midcutoffarray)):
    midarraydata.append([])
    
print "successfully created midarraydata"

# read the probeset file, which correlates name to position.
infilename = "probe_set_build37_forpartek.txt"
infile = open(infilename,"r")
indata = infile.readlines()
infile.close()

labels = {}
for line in indata:
    line = line.rstrip().split()
    if (len(line) == 3):
        (id, chr, pos) = line[0:3]
        if chr=="X":
            continue
#            chr = "23"
        if chr=="Y":
            continue
#            chr = "24"
        labels[id] = (chr, pos)

#Now read the SNP data file
samples = ["1047", "163", "222", "230", "312", "391", "611", "660", "729", "732", "824", "930"]
#samples = ["1047"]

for sample in samples:
    #Read the SNP data file
    SNPfilename = "SNP_data_orig/" + sample + "_BAFs"
    infile = open(SNPfilename,"r")
    indata = infile.readlines()
    infile.close()

    SNPs = {}
    labelrow = indata[0].rstrip().split("\t")
    BErow = indata[1].rstrip().split("\t")
    bloodrow = indata[2].rstrip().split("\t")
    for col in range(len(labelrow)):
        if labelrow[col] in labels:
            if (labelrow[col][0:3] == "cnv"):
                continue
            SNPs[labels[labelrow[col]]] = (BErow[col], bloodrow[col])
            if (BErow[col] != "?"):
                val = float(BErow[col])
                allarraydata.append(val)
                for c in range(len(midcutoffarray)):
                    if (val > midcutoffarray[c] and val < 1-midcutoffarray[c]):
                        midarraydata[c].append(val)
            if (bloodrow[col] != "?"):
                val = float(bloodrow[col])
                allarraydata.append(val)
                for c in range(len(midcutoffarray)):
                    if (val > midcutoffarray[c] and val < 1-midcutoffarray[c]):
                        midarraydata[c].append(val)
        

med = numpy.median(allarraydata)
print "median of all data = ", med
for c in range(len(midcutoffarray)):
    midmed = numpy.median(midarraydata[c])
    print "\nmedian of data greater than ", midcutoffarray[c], " from either 0 or 1: ", midmed
    
#OUTPUT, June 8th, 2016:
#median of all data =  0.53706
#median of data greater than  0.02  from either 0 or 1:  0.519243
#median of data greater than  0.04  from either 0 or 1:  0.522277
#median of data greater than  0.06  from either 0 or 1:  0.523133
#median of data greater than  0.08  from either 0 or 1:  0.523423
#median of data greater than  0.1  from either 0 or 1:  0.523576
#median of data greater than  0.12  from either 0 or 1:  0.523668
#median of data greater than  0.14  from either 0 or 1:  0.52372
#median of data greater than  0.16  from either 0 or 1:  0.523743
#median of data greater than  0.18  from either 0 or 1:  0.523747
#median of data greater than  0.2  from either 0 or 1:  0.523739
