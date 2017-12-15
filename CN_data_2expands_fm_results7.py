# -*- coding: utf-8 -*-
"""
Created on Wed May 11 17:10:20 2016

@author: lpsmith
"""

# convert kenji's CN calls (from exome data, IIRC) to actual CN values, in the format that EXPANDS wants.

comparefile = open("CN_compare.txt","w")
comparefile.write("chrom\tstartpos\tendpos\tCN_partek\tCN_kenji\n")

#read the kenji estimates
kenjifile = open("combined_results7.merged_calls", "r")
kenjidata = kenjifile.readlines()
ID_list = set()
outfile = kenjifile
for line in kenjidata:
    line = line.rstrip().split()
    if (line[0] == "#"):
        continue
    (id, chr, start, end, size, cn_call) = line[0:6]
    if (id not in ID_list):
        ID_list.add(id)
        outfile = open("data/" + id + "_CN_from_kenji_results7.txt", "w")
        outfile.write("chr\tstartpos\tendpos\tCN_Estimate\n")
    outfile.write(chr + "\t" + start + "\t" + end + "\t")
    cn_call = int(cn_call)
    if (cn_call == 0):
        outfile.write("0")
    elif(cn_call == 1):
        outfile.write("1")
    elif(cn_call == 2):
        outfile.write("2")
    elif(cn_call == 3):
        outfile.write("2")
    elif(cn_call == 4):
        outfile.write("3")
    elif(cn_call == 5):
        outfile.write("4")
    elif(cn_call == 6):
        outfile.write("6")
    else:
        outfile.write("UNKNOWN")
    outfile.write("\n")
