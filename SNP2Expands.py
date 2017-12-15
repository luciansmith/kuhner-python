# -*- coding: utf-8 -*-
"""
Created on Wed May 11 17:58:38 2016

@author: lpsmith
"""

# This file reads in SNP array BAF data, corrects for dye bias, and writes
# it out again in EXPANDS format.  There are basic filters that can be
# used, and it also can be filtered against the exome data.  Files with
# no filter, the basic filters, and basic+exome are all written out.


dye_median = 0.5237            #from running median_calc.py; June 8th, 2016
blood_be_diff_threshold = 0.03 #from running FilterFinder.py; June 17; 2016
het_threshold = 0.18           #from running FilterFinder_blood_call.py; June 17, 2016  
ex_snp_diff_threshold = 0.4    #from running FilterFinder_exome_diff.py; June 17, 2016
dbcorr = True


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

    #Now go through the SNP data and export different files for different filters
    dbcorrstr = ""
    if (dbcorr):
        dbcorrstr = "_dbcorr"
        
    outfn_0f = "SNP_data_filtered/" + sample + dbcorrstr + "_diff_filter.txt"
    outfn_1f = "SNP_data_filtered/" + sample + dbcorrstr + "_basic_filters.txt"
    outfn_2f = "SNP_data_filtered/" + sample + dbcorrstr + "_basic_and_exome_filters.txt"
    outf_0f = open(outfn_0f, "w")
    outf_1f = open(outfn_1f, "w")
    outf_2f = open(outfn_2f, "w")
    outline = "chr\tstartpos\tendpos\tAF_Tumor\tPN_B\tAF_Blood\n"
    outf_0f.write(outline)
    outf_1f.write(outline)
    outf_2f.write(outline)

    for label in SNPs.items():
        (chr, pos) = label[0]
        (SNP_BE, SNP_blood) = label[1]
        if (SNP_BE == "?" or SNP_blood == "?"):
            continue
        SNP_BE = float(SNP_BE)
        SNP_blood = float(SNP_blood)
        
        snpid = rev_labels[label[0]]
        if (snpid[0:3] != "cnv" and dbcorr):
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
                
        #Call for whether blood is heterozygous or not:
        ishet = "0"
        if (SNP_blood > 0.25 and SNP_blood < 0.75):
            ishet = "1"
        #If the two are not significantly different, don't write it out anywhere.
        if (abs(SNP_blood - SNP_BE) < blood_be_diff_threshold):
            continue
        
        #if the *called* blood frequency isn't significantly different from BE, don't write that out either.
        if (ishet == "0"):
            if (SNP_BE > 1-blood_be_diff_threshold or SNP_BE < blood_be_diff_threshold):
                continue
        else:
            if (abs(SNP_BE-0.5) < blood_be_diff_threshold):
                continue

        #If need be, flip which allele is 'B' so that the blood is the lower frequency
        if (SNP_blood > SNP_BE):
            SNP_BE = 1-SNP_BE
            SNP_blood = 1-SNP_blood

        #Create the output line that we may write out:
        outline = chr + "\t" + pos + "\t" + pos + "\t" + str(SNP_BE) + "\t" + ishet + "\t" + str(SNP_blood) + "\n"
        
        #Write this out to the no-filter file
        outf_0f.write(outline)
        

        #if the blood cannot be called as heterozygous or homozygous, don't write that out to a filtered file:
        if (SNP_blood > 0.0+het_threshold and SNP_blood < 0.5-het_threshold):
            continue
        if (SNP_blood > 0.5+het_threshold and SNP_blood < 1.0-het_threshold):
            continue
        
        #We're done with our first filter; write this out.
        outf_1f.write(outline)
        
        #now try to find the same data in the exome:
        exome = exomes.get(label[0])
        if (exome == None):
            continue
        
        (exomeBE, exomeBlood) = exome
        exomeBE = float(exomeBE)
        exomeBlood = float(exomeBlood)
        #if the two measures signifiantly differ, don't write them out:
        if (abs(SNP_BE - exomeBE) + abs(SNP_blood - exomeBlood) > ex_snp_diff_threshold):
            continue

        #This was our last filter; write it out.
        outf_2f.write(outline)
