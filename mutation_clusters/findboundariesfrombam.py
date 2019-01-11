#! /usr/bin/env python2

import sys
import pysam

min_mapping_quality = 25

def findminposin(readlist):
  minpos = None
  for read in readlist:
    if read.mapping_quality < min_mapping_quality:
      continue
    allpos = pysam.AlignedSegment.get_reference_positions(read,full_length=False)
    allpos.sort()
    if len(allpos) == 0:
      continue
    if minpos is None:
      minpos = allpos[0]
    if allpos[0] < minpos:
      minpos = allpos[0]

  return minpos

def findmaxposin(readlist):
  maxpos = None
  for read in readlist:
    if read.mapping_quality < min_mapping_quality:
      continue
    allpos = pysam.AlignedSegment.get_reference_positions(read,full_length=False)
    allpos.sort()
    if len(allpos) == 0:
      continue
    if maxpos is None:
      maxpos = allpos[-1]
    if allpos[-1] > maxpos: 
      maxpos = allpos[-1]

  return maxpos

def searchbackfrom(bamfile,chr,pos):
  oldminpos = -1
  minpos = pos
  while oldminpos != minpos:
    oldminpos = minpos
    reads = list(bamfile.fetch(chr,oldminpos,oldminpos+1))
    minpos = findminposin(reads)
    if minpos is None:
      return oldminpos

  return minpos

def searchforwardfrom(bamfile,chr,pos):
  oldmaxpos = -1
  maxpos = pos
  while oldmaxpos != maxpos:
    oldmaxpos = maxpos
    reads = list(bamfile.fetch(chr,oldmaxpos,oldmaxpos+1))
    maxpos = findmaxposin(reads)
    if maxpos is None:
      return oldmaxpos

  return maxpos

# search until a non-None value is found
def searchforwardby1k(bamfile,chr,pos):
  oldmaxpos = pos
  maxpos = None
  while maxpos is None:
    reads = list(bamfile.fetch(chr,oldmaxpos,oldmaxpos+1))
    maxpos = findmaxposin(reads)
    oldmaxpos += 1000

  return searchbackfrom(bamfile,chr,maxpos)

# search until a non-None value is found
def searchbackwardby1k(bamfile,chr,pos):
  oldminpos = pos
  minpos = None
  while minpos is None:
    reads = list(bamfile.fetch(chr,oldminpos,oldminpos+1))
    minpos = findminposin(reads)
    oldminpos -= 1000

  return searchforwardfrom(bamfile,chr,minpos)

# decide type and directionality of future searches
def searchin(bamfile,chr,pos,searchback):
  if pos == -1:
    return -1

  reads = list(bamfile.fetch(chr,pos,pos+1))
  if searchback:
    newpos = findminposin(reads)
  else:
    newpos = findmaxposin(reads)

  if newpos is None:
    # conduct a reverse search
    if searchback:
      return searchforwardby1k(bamfile,chr,pos)
    else:
      return searchbackwardby1k(bamfile,chr,pos)

  # conduct a normal search
  if searchback:
    return searchbackfrom(bamfile,chr,pos)
  else:
    return searchforwardfrom(bamfile,chr,pos)
  

############################################################################

if len(sys.argv) != 2:
  print "USAGE:  findboundariesfrombam.py bamfileurl"
  print "This version takes its BAM/BAI from the S3 cloud, and writes to ."
  exit()

bamurl = sys.argv[1]
items = bamurl.split("/")
pid,sid,dna,level = items[-1].split("-")

baiurl = bamurl[0:-1] + "i"
bamfile = pysam.AlignmentFile(bamurl,"rb", index_filename=baiurl)

startpos = {}
for line in open("chrom_boundaries","r"):
  line = line.rstrip()
  chr,startp,endp,startq,endq = line.split()
  startpos[chr] = [int(startp),int(endp),int(startq),int(endq)]

newpos = {}
for chr in startpos:
  allpos = startpos[chr]
  print "working on chromosome",chr
  minstartp = searchin(bamfile,chr,allpos[0],True)
  maxendp = searchin(bamfile,chr,allpos[1],False)
  minstartq = searchin(bamfile,chr,allpos[2],True)
  maxendq = searchin(bamfile,chr,allpos[3],False)
  newpos[int(chr)] = [minstartp,maxendp,minstartq,maxendq]
  print "finished chromesome",chr

delim = "_"
outfilename = pid+delim+sid+delim+dna+delim+level+delim
outfilename += "bam_chrom_boundaries.tsv"
outfile = open(outfilename,"w")
outline = "Patient\tSample\tChromosome\tpstart\tpend\tqstart\tqend\n"
outfile.write(outline)
for chr in sorted(newpos.iterkeys()):
  outline = pid+"\t"+sid+"\t"+str(chr)
  for pos in newpos[chr]:
    outline += "\t" + str(pos)
  outline += "\n"
  outfile.write(outline)
