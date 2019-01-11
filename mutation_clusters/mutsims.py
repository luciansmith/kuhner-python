# mutsims.py            Mary Kuhner and Jon Yamato 2018/04/26

# This program takes previously pickled information about mutation
# positions and uses it to build simulations, to which the real
# data from each sample can be comed.  Currently mutations are 2+
# caller NYGC Union files.

# Currently it does not filter out SCA regions.

import os, numpy
import random
import matplotlib.pyplot as plt
import pickle

nreps = 1000
limit = 50
allowed_chroms = [str(x) for x in range(1,23)]

########################################################################
# functions

def getdistances(muts):     # assumes muts is sorted!
  muts.sort()
  dists = []
  prevmut = muts[0]
  for mut in muts[1:]:
    dists.append(mut - prevmut)
    prevmut = mut
  return dists

def bincount(muts,limit):
  dists = getdistances(muts)
  result = [0,0]
  for dist in dists:
    if dist <= limit:
      result[0] += 1
    else:
      result[1] += 1
  return result

def simulate_data(mutstore,target,limit):
# randomly choose "target" mutations from the real data for this arm
  assert len(mutstore) >= target    # not enough mutations, uh-oh!
  allbins = []
  for rep in range(nreps):
    mutset = set()
    while len(mutset) < target:            # we don't want duplicates
      mutset.add(random.choice(mutstore))
    muts = list(mutset)
    muts.sort()
    bins = bincount(muts,limit)
    allbins.append(bins)
  return allbins

def indexscore(rbins,sbins):
  r_p50 = float(rbins[0])/sum(rbins)
  s_p50 = [float(x[0])/sum(x) for x in sbins]
  s_p50.sort()
  tie_indexes = []
  for i in range(nreps):
    if r_p50 > s_p50[i]:
      continue
    if r_p50 == s_p50[i]:
      tie_indexes.append(i)
      continue
    if r_p50 < s_p50[i]:  
      # if we never found a tie, the previous value was the answer
      # otherwise the ties are the answer
      if len(tie_indexes) == 0:  
        tie_indexes.append(i-1)
      break
  numties = len(tie_indexes)
  if numties == 0:       # real value bigger than all simulated
    return nreps
  elif numties == 1:  # a single answer
    return tie_indexes[0]
  else:
    target = numties/2    # ties -- yes, integer math
    return tie_indexes[target]

def accumulate_overall(rbins,rbins_overall):
  for i in range(len(rbins)):
    rbins_overall[i] += rbins[i]

def accumulate_overall_sim(sbins,sbins_overall):
  for rep in range(nreps):
    accumulate_overall(sbins[rep],sbins_overall[rep])

def simscore(nreps,arms,muts,mutstore,limit,arm_scores,overall_scores):
   # determine distances
   rbins_overall = [0,0]
   sbins_overall = [[0,0] for x in range(nreps)]
   for arm in arms:
     if arm not in muts:    # no mutations on this arm in this sample
       continue
     nmuts = len(muts[arm])
     if nmuts < 2:   # not enough mutations on this arm in this sample
       continue

     # calculate for real data
     rbins = bincount(muts[arm],limit)
     # simulate
     sbins = simulate_data(mutstore[arm],nmuts,limit)
     # calculate score for this arm
     score = indexscore(rbins,sbins)
     arm_scores[arm].append(score)
     # accumulate values to calculate score over whole genome
     accumulate_overall(rbins,rbins_overall)
     accumulate_overall_sim(sbins,sbins_overall)

   if sum(rbins_overall) == 0:
     print("No mutations in this sample!?")
   else:
     overall_scores.append(indexscore(rbins_overall,sbins_overall))


########################################################################
# main program

# read chromosome boundaries file
boundfile = "/home/mkkuhner/seagate/data/wgs_dna/bam/chrom_boundaries"

boundaries = {}
arm_boundaries = {}
#filepath = boundlocation + file
for line in open(boundfile,"r"):
  chr,pstart,pend,qstart,qend = line.rstrip().split()
  pstart = int(pstart)
  pend = int(pend)
  qstart = int(qstart)
  qend = int(qend)
  boundaries[chr] = [pstart,pend,qstart,qend]
  if pstart != -1:
    armname = chr+"p"
    arm_boundaries[armname] = [pstart,pend]
  armname = chr+"q"
  arm_boundaries[armname] = [qstart,qend]

arms = arm_boundaries.keys()


# unpickle stored mutation location data
pfile = open("mutpositions.pkl","r")
mutstore_overall = pickle.load(pfile)
mutstore_mod = pickle.load(pfile)
mutstore_nmod = pickle.load(pfile)


# locate mutation files
mutlocation = "/home/mkkuhner/seagate/data/wgs_dna/union/"
mutfiles = []
for root, dirs, files in os.walk(mutlocation):
  for file in files:
    if file.endswith(".annotated.txt"):
      mutfiles.append(file)
print("Assessing mutations in",len(mutfiles),"input files")

# initialize variables
overall_scores = []
arm_scores = {}
for arm in arms:
  arm_scores[arm] = []

overall_scores_mod = []
arm_scores_mod = {}
for arm in arms:
  arm_scores_mod[arm] = []

overall_scores_nmod = []
arm_scores_nmod = {}
for arm in arms:
  arm_scores_nmod[arm] = []

# process mutation files
for mutfile in mutfiles:   
   muts = {}
   muts_mod = {}
   muts_nmod = {}
   prefix = mutfile.split("/")[-1]
   prefix = prefix.split("-")
   pid = prefix[0]
   sid = prefix[1]
   prefix = pid + "-" + sid + "-"
   samplename = pid+"_"+sid

   if sid.endswith("N"):
     print("Skipping",pid,sid)
     continue
  
   print("Analyzing",pid,sid)
 
   #################################
   # parse the mutation file to find eligible mutation pairs
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
     call = call.split("-")
     # if call != "mutect-lofreq-strelka_snv":  # require all three callers
     #   continue
     if len(call) < 2:   # require at least two callers
       continue
     lstart,lend,rstart,rend = boundaries[chr]
     armname = None
     if pos >= lstart and pos <= lend:
       armname = chr+"p"
     elif pos >= rstart and pos <= rend:
       armname = chr+"q"
     if armname == None:
       # print("stray mutation:",chr,lstart,lend,rstart,rend, "don't cover",pos)
       continue 
     else:
       count += 1

     if armname not in muts:
       muts[armname] = []
     muts[armname].append(pos)

     snpeffcall = line[snpeffindex]
     snpeffcalls = ["HIGH","MODERATE","LOW","MODIFIER"]
     if snpeffcall not in snpeffcalls:
       print("\nfound an abberant snpeffcall for chomosome",chr,)
       print("position",pos,":  ",snpeffcall)
       exit()
     if snpeffcall == snpeffcalls[3]:
       if armname not in muts_mod:
         muts_mod[armname] = []
       muts_mod[armname].append(pos)
     else: # was something other than "modifier"
       if armname not in muts_nmod:
         muts_nmod[armname] = []
       muts_nmod[armname].append(pos)

   simscore(nreps,arms,muts,mutstore_overall,limit,arm_scores,overall_scores)
   print("overall",)
   simscore(nreps,arms,muts_mod,mutstore_mod,limit,arm_scores_mod,overall_scores_mod)
   print("modifier",)
   simscore(nreps,arms,muts_nmod,mutstore_nmod,limit,arm_scores_nmod,overall_scores_nmod)
   print("non-modifier")

pfile = open("indexscores.pkl","w")
pickle.dump(overall_scores,pfile)
pickle.dump(arm_scores,pfile)
pickle.dump(overall_scores_mod,pfile)
pickle.dump(arm_scores_mod,pfile)
pickle.dump(overall_scores_nmod,pfile)
pickle.dump(arm_scores_nmod,pfile)
pfile.close()

#import matplotlib.pyplot as plt

#figno = 0
#plt.figure(figno)
#plt.title("Overall results")
#plt.hist(overall_scores)

#for arm in arms:
#  figno += 1
#  plt.figure(figno)
#  plt.title("Chrom"+arm)
#  plt.hist(arm_scores[arm])
  
#plt.show()
