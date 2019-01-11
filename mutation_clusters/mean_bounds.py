import os

indir = "results/"
infiles = []
for root,dirs,files in os.walk("results/"):
  for file in files:
    if file.endswith("tsv"):
      infiles.append(file)

chromdict = {}
numfiles = 0.0
for infile in infiles:
  numfiles += 1.0
  for line in open(indir + infile,"r"):
    if line.startswith("Patient"):    # skip header
      continue
    line = line.rstrip().split()
    chrom = line[2]
    pvals = line[3:7]
    pvals = [int(x) for x in pvals]
    if chrom not in chromdict:
      chromdict[chrom] = [0.0 for x in range(4)]
    for i in range(4):
      if pvals[i] != -1:
        chromdict[chrom][i] += pvals[i]
      else:
        chromdict[chrom][i] = -1

for chrom in chromdict.keys():
  for i in range(len(chromdict[chrom])):
    if chromdict[chrom][i] != -1:
      chromdict[chrom][i] /= numfiles

sortkeys = chromdict.keys()
sortkeys = [int(x) for x in sortkeys]
sortkeys.sort()
sortkeys = [str(x) for x in sortkeys]

outfile = open("new_chrom_boundaries","w")
outline = "Patient\tSample\tChromosome\tpstart\tpend\tqstart\tqend\n"
outfile.write(outline)
for i in range(len(sortkeys)):
  chrom = sortkeys[i]
  outline = chrom 
  for j in range(4):
    outline += "\t" + str(int(chromdict[chrom][j])) 
  outline += "\n"
  outfile.write(outline) 
