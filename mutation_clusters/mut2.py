# mutdistribution.py            Mary Kuhner and Jon Yamato 2018/04/26

# This program finds the distribution of distances among mutations in a 
# mutation file (currently the Union files from NYGC, but could easily
# work with a strelka .vcf instead).

# Currently it does not filter out SCA regions.

import os, numpy
import random


########################################################################
# functions

def makerandommutations(nmut,posmax):
  mutsat = random.sample(xrange(posmax),nmut)
  mutsat.sort()
  return mutsat

def makedistancepairs(muts):
  distpairs = []
  prevmut = muts[0]
  for mut in muts[1:]:
    distpairs.append((prevmut,mut))
    prevmut = mut
  return distpairs    

def simulatemutations(actualstd,target,nmut):
  nreps = 100
   
  #print("sprinkling",nmut,"mutations on",mbtarget,"megabases",nreps,"times")

  stds = []
  for rep in xrange(nreps):
    mutsat = makerandommutations(nmut,target)
    distpairs = makedistancepairs(mutsat)
  
    distancebetween = []
    for m1,m2 in distpairs:
      distancebetween.append(m2-m1)
  
    newstd = numpy.std(distancebetween,ddof=1)
  
    stds.append(newstd)

  stds.sort()
  for i in xrange(len(stds)):
    if actualstd < stds[i]:
      return float(i)/nreps
    

########################################################################
# main program

allowed_chroms = [str(x) for x in range(1,23)]

# read chromosome boundaries file
boundaries = {}
arm_boundaries = {}
for line in open("chrom_boundaries","r"):
  chr,lstart,lend,rstart,rend = line.rstrip().split()
  lstart = int(lstart)
  lend = int(lend)
  rstart = int(rstart)
  rend = int(rend)
  boundaries[chr] = [lstart,lend,rstart,rend]
  if lstart != -1:
    armname = chr+"p"
    arm_boundaries[armname] = [lstart,lend]
  armname = chr+"q"
  arm_boundaries[armname] = [rstart,rend]


# process mutation files
mutlocation = "/home/mkkuhner/seagate/data/wgs_dna/union/"
mutfiles = []
for root, dirs, files in os.walk(mutlocation):
  for file in files:
    if file.endswith(".annotated.txt"):
      mutfiles.append(file)

sigs = {}

for mutfile in mutfiles:   
   muts = {}
   prefix = mutfile.split("/")[-1]
   prefix = prefix.split("-")
   pid = prefix[0]
   sid = prefix[1]
   prefix = pid + "-" + sid + "-"
   print("Analyzing",pid,sid)
 
   #################################
   # parse the mutation file to find eligible mutation pairs
   prevch = 0
   prevpos = -500
   
   first = True 
   count = 0
   for line in open(mutlocation + mutfile,"r"):
     line = line.rstrip().split()
     # process header
     if first == True:
       chromindex = line.index("CHROM")
       posindex = line.index("POS")
       refindex = line.index("REF")
       altindex = line.index("ALT")
       callindex = line.index("CALLED_BY")
       snpeffindex = line.index("SNPEFF_IMPACT")
       first = False
       continue    

     # process mutation
     chr = line[chromindex]
     if chr not in allowed_chroms:   # no X, Y, mt, etc.
       continue
     pos = int(line[posindex])
     call = line[callindex]
     if call != "mutect-lofreq-strelka_snv":  # require all three callers
       continue
     lstart,lend,rstart,rend = boundaries[chr]
     armname = None
     if pos >= lstart and pos <= lend:
       armname = chr+"p"
     if pos >= rstart and pos <= rend:
       armname = chr+"q"
     if armname == None:
       # print("stray mutation:",chr,lstart,lend,rstart,rend, "don't cover",pos)
       continue 
     else:
       count += 1

     if armname not in muts:
       muts[armname] = []
     muts[armname].append(pos)

   # determine distances
   distances = {}
   for arm in muts.keys():
     distances[arm] = []
     prevpos = -1
     for pos in muts[arm]:
       if prevpos == -1:
         prevpos = pos
         continue
       distances[arm].append(pos - prevpos)
       prevpos = pos

   # assess significance for each arm
   for arm in distances.keys():
     if len(muts[arm]) < 4:
       continue
     astart,aend = arm_boundaries[arm]
     alength = aend - astart + 1
     sortmuts = muts[arm]
     if len(sortmuts) <= 1:  continue
     nummuts = len(sortmuts)
     sortmuts.sort()
     actualstd = numpy.std(sortmuts,ddof=1)
     sig = simulatemutations(actualstd,alength,nummuts)
     if arm not in sigs:
       sigs[arm] = []
     sigs[arm].append(sig)

print("Arm\tTotal\tSignificant")
for arm in sigs.keys():
  # count how many sigs are >= 0.95
  count = 0
  for sig in sigs[arm]:
    if sig >= 0.95:
      count += 1
  print(arm, len(sigs[arm]),count)
