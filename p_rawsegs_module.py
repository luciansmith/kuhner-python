# take the "raw segmentation" dump from P-ASCAT and APPEND it to an
# existing all_lesions.txt file (create file if it doesn't exist).

# relies on Xiaohong's list of chromosome boundaries
# found in "chrom_boundaries"

# CONVENTIONS:  Throughout this program "chr" the chromosome ID is
# always a STRING, while pos, startpos, endpos etc. are INTS.

# We pull copy numbers from the ID_raw_segments.txt file from 
# P-ASCAT and cross-reference the starting and ending probes with
# the probeset to determine segment boundaries and chromosomes.

validchroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
  "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]

import numpy

# DEBUG function to print segs of a chrom
def printsegs(chr, segsbychrom):
  for seg in segsbychrom[chr]:
    print seg

# find the next SNP in the probeset, or None if there isn't one
def nextsnp(chr,pos,logRbychrom):
  data = logRbychrom[chr]
  for i in xrange(len(data)):
    if data[i][1] == pos:
      if i+1 != len(data):
        return data[i+1][1]
      else:
        return None

# find the previous SNP in the probeset, or None if there isn't one
def prevsnp(chr,pos,logRbychrom):
  data = logRbychrom[chr]
  for i in xrange(len(data)):
    if data[i][1] == pos:
      if i-1 >= 0:
        return data[i+1][1]
      else:
        print "prestart SNP at",chr,pos
        return None

# remove # and _ from SNP names before use; R does not like them!
def correctname(name):
  if "#" in name:
    name = name.split("#")
    newname = ""
    for part in name:
      newname += part
    name = newname
  if "_" in name:
    name = name.split("_")
    newname = ""
    for part in name:
      newname += part
    name = newname
  return name

def readblock(data, keywd):
  # find chromosomes
  collection = []
  found = False
  done = False
  for line in data:
    if done: break     # we're past desired section

    line = line.rstrip()
    line = line.rstrip(",")
    line = line.replace(" ","")
    line = line.replace('"','')

    # the line consists of quoted comma-delimited entries with spaces
    # between them; a whole section is enclosed in parentheses.  

    if keywd in line:  # start of section
      line = line.split("(")  
      line = line[-1]         
      found = True
      line = line.rstrip().split(",")   
      for entry in line:
        collection.append(entry)
      found = True
      continue

    if found and ")" not in line:   # middle of section
      line = line.split(",")
      for entry in line:
        collection.append(entry)
      continue

    if found:
      line = line.split(")")
      line = line[0]
      if line != "":     # test for the lonely parenthesis
        line = line.split(',')
        for entry in line:
          collection.append(entry)
      done = True

  return collection

def makeprobedict():
  # path to file that explains chip layout
  probesetname = "PATH_TO_DATA/probe_set_build37_forpartek.txt"
  # set up dictionary relating probe name to chromosome and position
  probes = {}
  first = True
  for line in open(probesetname,"r"):
    if first:   # skip header
      first = False
      continue
    line = line.rstrip()
    line = line.split()
    mykey = line[0]
    mykey = correctname(mykey)
    mydata = (line[1],int(line[2]))
    probes[mykey] = mydata
  return probes

def getlogR(chrom,startpos,endpos,logRbychrom):
  # find the logR for the segment at (chrom,startpos-endpos)
  # by averaging SNP logR values
  logRvals = []
  for snp in logRbychrom[chrom]:
    if snp[1] >= startpos and snp[1] <= endpos:
      logRvals.append(snp[2])
  return numpy.mean(logRvals)

def distbetween(chrom,startpos,endpos,snpsbychrom):
  starti = snpsbychrom[chrom].index(startpos)
  endi = snpsbychrom[chrom].index(endpos)
  return endi - starti + 1

class lesion:
  def __init__(self,chr,startpos,endpos,rawA,rawB,intA,intB,logR):
    self.chr = str(chr)
    self.startpos = int(startpos)
    self.endpos = int(endpos)
    self.rawA = rawA
    self.rawB = rawB
    if rawB != None:
      self.rawcn = float(rawA) + float(rawB)
    else:
      self.rawcn = None
    if intB != None:
      self.intA = intA
      self.intB = intB
      self.intcn = int(intA) + int(intB)
    else:
      self.intcn = intA
      self.intA = self.intcn/2
      self.intB = self.intcn/2
    if logR != None:
      self.logR = float(logR)
    else:
      self.logR = None
    self.index = None
  def __str__(self):
    outline = self.chr + " " + str(self.startpos) + " " + str(self.endpos) 
    outline += " int CN " +  str(self.intcn) + " raw CN " + str(self.rawcn)
    outline += " logR " + str(self.logR) + " index " + str(self.index)
    return outline

def chrom_start_end_mapper(directory,id):
  copyfilename = directory+"/"+id + "_copynumber_segments.txt"
  copyfile = open(copyfilename,"r")
  copyfile.readline()    # discard header
  cmap = {}
  for line in copyfile.readlines():
    line = line.split()
    chr,st,en = line[0:3]
    if (st,en) not in cmap:
      print "Mapping",st,en,"to chrom",chr
      cmap[(st,en)] = chr
    else:
      # is this just the same exact segment again?
      if chr == cmap[(st,en)]:  continue
      # Uh-oh, it is NOT
      print "Non-unique start/end pair detected"
      print "Ends are",st,en
      print "Found on, at least,",cmap[(st,en)],chr
      exit()
  return cmap

##### main program #####

