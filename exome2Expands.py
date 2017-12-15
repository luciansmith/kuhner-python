# convert the raw exome sequencing BAF data to EXPANDS-ready input

from os import walk
from os import path

blood_be_diff_threshold = 0.03 # From ExomeFilterFinder_bloodBEdiff.py
het_threshold = 0.18           # From ExomeFilterFinder_bloodcall.py
ex_snp_diff_threshold = 0.4    # From FilterFinder_exome_diff.py (same as for SNP)


# read the probeset file, which correlates name to position.
infilename = "probe_set_build37_forpartek.txt"
infile = open(infilename,"r")
indata = infile.readlines()
infile.close()

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

#filenames = ["exome_data_orig/1047_18966_180B_28_muTect.baf"]
filenames = []
for (_, _, f) in walk("exome_data_orig/"):
    filenames += f
    break        

for file in filenames:
    patientID = file[0:file.find("_")]
    
    #Read the SNP data file
    SNPs = {}
    SNPfilename = "SNP_data_orig/" + patientID + "_BAFs"
    hasSNPdata = False
    if (path.isfile(SNPfilename)):
        hasSNPdata = True
        infile = open(SNPfilename,"r")
        indata = infile.readlines()
        infile.close()
    
        #Process the SNP data
        labelrow = indata[0].rstrip().split("\t")
        BErow = indata[1].rstrip().split("\t")
        bloodrow = indata[2].rstrip().split("\t")
        for col in range(len(labelrow)):
            if labelrow[col] in labels:
                SNPs[labels[labelrow[col]]] = (BErow[col], bloodrow[col])
                #print labels[labelrow[col]]

    # read a Patty logr file
    infile = open("exome_data_orig/" + file,"r")
    indata = infile.readlines()
    infile.close()

    outfn_0f = "exome_data_filtered/" + patientID + "_diff_filter.txt"
    outfn_1f = "exome_data_filtered/" + patientID + "_basic_filters.txt"
    outline = "chr\tstartpos\tendpos\tAF_Tumor\tPN_B\tAF_Blood\n"
    outf_0f = open(outfn_0f, "w")
    outf_1f = open(outfn_1f, "w")
    outf_0f.write(outline)
    outf_1f.write(outline)
    if (hasSNPdata):
        outfn_2f = "exome_data_filtered/" + patientID + "_basic_and_exome_filters.txt"
        outf_2f = open(outfn_2f, "w")
        outf_2f.write(outline)
 
    numfail = 0
    for lineorig in indata:
        line = lineorig.rstrip().split()
        try:
            (chr,name,pos,exomeBlood,exomeBE) = line[0:5]
        except ValueError:
            print "problematic line: from patient data " + patientID + ": '" + lineorig + "'"
            continue
        chr = chr[3:len(chr)]
        if (chr=="X"):
            chr = "23"
        if (chr=="Y"):
            chr = "24"
        if (chr[0]<'0' or chr[0]>'9'):
            continue
        blood_ishet = "0"
        try:
            exomeBlood = float(exomeBlood)
            exomeBE = float(exomeBE)
        except ValueError:
            print "Well, that was weird: " + lineorig
            continue
        if(exomeBlood > 0.25 and exomeBlood < 0.75):
            blood_ishet = "1"
        
        #If need be, flip which allele is 'B' so that the blood is the lower frequency
        if (exomeBlood > exomeBE):
            exomeBE = 1-exomeBE
            exomeBlood = 1-exomeBlood

        #Filter 1:  ensure that blood and BE are different from each other.
        if (abs(exomeBlood-exomeBE) > blood_be_diff_threshold):
            continue

        #if the *called* blood frequency isn't significantly different from BE, don't write that out either.
        if (blood_ishet == "0"):
            if (exomeBE > 1-blood_be_diff_threshold or exomeBE < blood_be_diff_threshold):
                continue
        else:
            if (abs(exomeBE-0.5) < blood_be_diff_threshold):
                continue

        #Create the output line that we may write out:
        outline = chr + "\t" + pos + "\t" + pos + "\t" + str(exomeBE) + "\t" + blood_ishet + "\t" + str(exomeBlood) + "\n"

        #That's the only filter we need for our first file:
        outf_0f.write(outline)
        
        #Filter 2:  check if the blood cannot be called as heterozygous or homozygous.
        if (exomeBlood > 0.0+het_threshold and exomeBlood < 0.5-het_threshold):
            continue
        if (exomeBlood > 0.5+het_threshold and exomeBlood < 1.0-het_threshold):
            continue
        
        #We're done with our 'called' filter; write this out.
        outf_1f.write(outline)
        
        if (hasSNPdata and rev_labels.get((chr, pos)) != None):
            id = rev_labels[(chr, pos)]
            SNP = SNPs.get(id)
            if (SNP == None):
                continue
            (SNPbe, SNPblood) = SNP
            SNPbe = float(SNPbe)
            SNPblood = float(SNPblood)
            #if the two measures signifiantly differ, don't write them out:
            if (abs(SNPbe - exomeBE) + abs(SNPblood - exomeBlood) > ex_snp_diff_threshold):
                continue

            #This was our last filter; write it out.
            outf_2f.write(outline)
            
