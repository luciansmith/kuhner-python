# compute unnormalized mutation signatures from trinucleotide context
# data, 2+ callers, from NYGC

# NB:  This uses the REF trinucleotide context, not the one from the patient's
# normal sample, alas.

# NB:  We assume that the input mutations are sorted by chromosome and position
from __future__ import division
from os import walk
from os import path
from os import mkdir
from random import shuffle

import lucianSNPLibrary as lsl

trinucdir = "trinucleotide/"
outdir = "tri_out/"


if not(path.isdir(outdir)):
    mkdir(outdir)

infiles = []
for root,dirs,files in walk(trinucdir):
  for file in files:
    if file.endswith("trinuc.txt"):
      infiles.append(file)

outfile = open(outdir + "mutfile_all.tsv","w")
outline = "sample.id\tchr\tpos\tref\talt\n"
outfile.write(outline)

mutcountfile = open(outdir + "mutcount_all.tsv", "w")
mutcountfile.write("Patient\tSample Group\tMutation count\n")

def getIdFrom(patient, samplelist):
    id = patient
    for sample in samplelist:
        id += "_" + sample
    return id

mutations = {}

# process mutation files
for file in infiles:
  patient,sample = file.split("-")[0:2]
  print("Processing",patient, ",", sample)
  if patient not in mutations:
      mutations[patient] = {}

  for line in open(trinucdir + file,"r"):
    line = line.rstrip().split()
    # process header line
    if line[0].startswith("CHROM"):  
      chromind = line.index("CHROM")
      posind = line.index("POS")
      trinucind = line.index("TRINUC")
      refind = line.index("REF")
      altind = line.index("ALT")
      continue

    # process non header line
    chr = line[chromind]
    pos = line[posind]
    alt = line[altind]
    ref = line[refind]
    mutation = (chr, pos, ref, alt)
    if mutation not in mutations[patient]:
        mutations[patient][mutation] = []
    mutations[patient][mutation].append(sample)

sorted_muts = {}
for patient in mutations:
    for mutation in mutations[patient]:
        newid = getIdFrom(patient, mutations[patient][mutation])
        if newid not in sorted_muts:
            sorted_muts[newid] = []
        sorted_muts[newid].append(mutation)
#        outfile.write(getIdFrom(patient, mutations[patient][mutation]))
#        outfile.write("\t" + mutation[0])
#        outfile.write("\t" + mutation[1])
#        outfile.write("\t" + mutation[2])
#        outfile.write("\t" + mutation[3])
#        outfile.write("\n")
#        if len(mutations[patient][mutation]) == 1:
#            outfile.write(patient + "_tip")
#            outfile.write("\t" + mutation[0])
#            outfile.write("\t" + mutation[1])
#            outfile.write("\t" + mutation[2])
#            outfile.write("\t" + mutation[3])
#            outfile.write("\n")

for newid in sorted_muts:
    countid = newid.replace("_","\t", 1)
    mutcountfile.write(countid + "\t" + str(len(sorted_muts[newid])) + "\n")
    if len(sorted_muts[newid]) < 50:
        continue
    for mutation in sorted_muts[newid]:
        outfile.write(newid)
        outfile.write("\t" + mutation[0])
        outfile.write("\t" + mutation[1])
        outfile.write("\t" + mutation[2])
        outfile.write("\t" + mutation[3])
        outfile.write("\n")
    
    

outfile.close()
