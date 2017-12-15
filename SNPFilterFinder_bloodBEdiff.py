# -*- coding: utf-8 -*-
"""
Created on Wed May 11 17:58:38 2016

@author: lpsmith
"""

# Based on SNP2Expands, but written to help find out what values that program
# should use for its filters.  Collects data for which we have both exome and
# SNP data, and writes out *all* of it.

dye_median = 0.5237 #from running median_calc.py; June 8th, 2016

# read the probeset file, which correlates name to position.
infilename = "probe_set_build37_forpartek.txt"
infile = open(infilename,"r")
indata = infile.readlines()
infile.close()

samples = ["1047", "163", "222", "230", "312", "391", "611", "660", "729", "732", "824", "930"]
#samples = ["1047"]

labels = {}
rev_labels = {}
for line in indata:
    line = line.rstrip().split()
    if (len(line) == 3):
        (id, chr, pos) = line[0:3]
        if chr=="X":
            chr = "23"
        if chr=="Y":
            chr = "24"
        labels[id] = (chr, pos)
        rev_labels[(chr, pos)] = id
        
outfilename = "filterfinder_SNP_bloodBEdiff.txt"
outfile = open(outfilename, "w")
outline = "Patient\tSNP_id\tchr\tpos\tSNP_BE\tSNP_blood\texomeBE\texomeBlood\terror\t0\t0.01\t0.02\t0.03\t0.04\t0.05\t0.06\t0.07\t0.08\t0.09\t0.1\t0.11\t0.12\t0.13\t0.14\t0.15\t0.16\t0.17\t0.18\t0.19\t0.2\n"
outfile.write(outline)

for sample in samples:
    #Read the SNP data file
    SNPfilename = "SNP_data_orig/" + sample + "_BAFs"
    infile = open(SNPfilename,"r")
    indata = infile.readlines()
    infile.close()
    
    #Process the SNP data
    SNPs = {}
    labelrow = indata[0].rstrip().split("\t")
    BErow = indata[1].rstrip().split("\t")
    bloodrow = indata[2].rstrip().split("\t")
    for col in range(len(labelrow)):
        if labelrow[col] in labels:
            SNPs[labels[labelrow[col]]] = (BErow[col], bloodrow[col])
            #print labels[labelrow[col]]

    #Process the exome data.
    exomes = {}
    exomefilename = "exome_data_reformatted/" + sample + "_full.baf"
    infile = open(exomefilename, "r")
    indata = infile.readlines()
    infile.close()
    for exline in indata:
        exline = exline.rstrip().split("\t")
        if exline[1] == "startpos":
            continue
        labelpair = (exline[0], exline[1])
        exomes[labelpair] = (exline[3], exline[5])

    for label in SNPs.items():
        (chr, pos) = label[0]
        (SNP_BE, SNP_blood) = label[1]
        if (SNP_BE == "?" or SNP_blood == "?"):
            continue
        SNP_BE = float(SNP_BE)
        SNP_blood = float(SNP_blood)
        
        snpid = rev_labels[label[0]]
        if (snpid[0:3] == "cnv"):
            continue
        #dye bias correction, BE:
        if (SNP_BE < dye_median):
            SNP_BE = SNP_BE * 0.5/dye_median
        else:
            SNP_BE = 0.5 + 0.5*(SNP_BE-dye_median)/(1-dye_median)

        #dye bias correction, blood:
        if (SNP_blood < dye_median):
            SNP_blood = SNP_blood * 0.5/dye_median
        else:
            SNP_blood = 0.5 + 0.5*(SNP_blood-dye_median)/(1-dye_median)
                
        #now try to find the same data in the exome:
        exome = exomes.get(label[0])
        if (exome == None):
            continue
        
        (exomeBE, exomeBlood) = exome
        exomeBE = float(exomeBE)
        exomeBlood = float(exomeBlood)
        
        #If need be, flip which allele is 'B' so that the blood is the lower frequency
        if (SNP_blood > SNP_BE):
            SNP_BE = 1-SNP_BE
            SNP_blood = 1-SNP_blood

        if (exomeBlood > exomeBE):
            exomeBE = 1-exomeBE
            exomeBlood = 1-exomeBlood
        #Now calculate the error, and whether it's filtered or not, and when.
        error = abs(exomeBlood-SNP_blood) + abs(exomeBE-SNP_BE)
        errorline = str(error)
        for fval in [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2]:
            if (abs(SNP_blood - SNP_BE) > fval):
                errorline += "\t" + str(error)

        #Create the output line
        outline = str(sample) + "\t" + snpid + "\t" + chr + "\t" + pos + "\t" + str(SNP_BE) + "\t" + str(SNP_blood) + "\t" + str(exomeBE) + "\t" + str(exomeBlood) + "\t" + errorline + "\n"
        outfile.write(outline)

print "done"