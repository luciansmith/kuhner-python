# mutdistribution.py            Mary Kuhner and Jon Yamato 2018/04/26

# This program finds the distribution of distances among mutations in a 
# mutation file (currently the Union files from NYGC, but could easily
# work with a strelka .vcf instead).

# Currently it does not filter out SCA regions.

import os, numpy
import random
import matplotlib.pyplot as plt

categories = 2000
nreps = 1000

########################################################################
# functions

def makerandommutations(nmut,posmax):
  mutsat = random.sample(xrange(posmax),nmut)
  mutsat.sort()
  return mutsat

def getdistances(muts):
  dists = []
  prevmut = muts[0]
  for mut in muts[1:]:
    dists.append(mut - prevmut)
  return dists

def categorize_distances(accumulator,distances,categories):
  # note:  accumulator must have categories+1 cells!
  assert len(accumulator) == categories + 1
  for dist in distances:
    if dist < categories:
      accumulator[dist] += 1.0
    else:
      accumulator[categories] += 1.0

def simulate_distances(target,nmut,arm):
  dist_overall = []
  dist_rep = [[] for x in xrange(repcount)]
  repcount = 0
  for rep in xrange(nreps):
    # note:  these mutations are pre-sorted 
    mutsat = makerandommutations(nmut,target)
    dist = getdistances(mutsat)
    dist_overall += dist
    dist_rep[repcount] = dist
    repcount += 1
  return dist_overall, dist_rep

def calc_distances(muts,arm_bounds):
  distances = {}
  r_distances = [0.0 for x in xrange(categories+1)]
  s_distances = [0.0 for x in xrange(categories+1)]
  for arm in muts.keys():
    # I assume these mutations are already sorted
    distances[arm] = getdistances(muts[arm])
    categorize_distances(r_distances,distances[arm],categories)

    # simulated distances
    astart,aend = arm_bounds[arm]
    alength = aend - astart + 1
    nummuts = len(muts[arm])
    if nummuts < 2:  continue     # can't make distances!
    sim_overall, sim_rep = simulate_distances(alength,nummuts,arm)
    categorize_distances(s_distances,sim_overall,categories)

  return r_distances,s_distances

def calcperc(countlist,limit):
  total = float(sum(countlist))
  if total == 0.0:
    return 0.0
  #using "limit+1" because we want to include the limit'th entry.
  return sum(countlist[:limit+1]) / total

########################################################################
# main program

allowed_chroms = [str(x) for x in range(1,23)]

# read chromosome boundaries file
boundfile = "/home/mkkuhner/seagate/data/wgs_dna/bam/chrom_boundaries"

#boundfiles = []
#for root, dirs, files in os.walk(boundlocation):
#  for file in files:
#    if file.endswith("_boundaries.tsv"):
#      boundfiles.append(file)
#print "Reading boundaries from",len(boundfiles),"files"

#sample_boundaries = {}
#for file in boundfiles:

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

# process mutation files
mutlocation = "/home/mkkuhner/seagate/data/wgs_dna/union/"
mutfiles = []
for root, dirs, files in os.walk(mutlocation):
  for file in files:
    if file.endswith(".annotated.txt"):
      mutfiles.append(file)
print "Assessing mutations in",len(mutfiles),"input files"

sigs = {}
sim_distances = [0.0 for x in xrange(categories+1)]
real_distances = [0.0 for x in xrange(categories+1)]
mod_sim_distances = [0.0 for x in xrange(categories+1)]
mod_real_distances = [0.0 for x in xrange(categories+1)]
nmod_sim_distances = [0.0 for x in xrange(categories+1)]
nmod_real_distances = [0.0 for x in xrange(categories+1)]
real_50count = []
sim_50count = []

for mutfile in mutfiles:   
   muts = {}
   modifiercalls = {}
   notmodifiercalls = {}
   prefix = mutfile.split("/")[-1]
   prefix = prefix.split("-")
   pid = prefix[0]
   sid = prefix[1]
   prefix = pid + "-" + sid + "-"
   samplename = pid+"_"+sid

   if sid.endswith("N"):
     print "Skipping",pid,sid
     continue
  
   print "Analyzing",pid,sid
 
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
       # print "stray mutation:",chr,lstart,lend,rstart,rend, "don't cover",pos
       continue 
     else:
       count += 1

     if armname not in muts:
       muts[armname] = []
     muts[armname].append(pos)

     snpeffcall = line[snpeffindex]
     snpeffcalls = ["HIGH","MODERATE","LOW","MODIFIER"]
     if snpeffcall not in snpeffcalls:
       print "\nfound an abberant snpeffcall for chomosome",chr,
       print "position",pos,":  ",snpeffcall
       exit()
     if snpeffcall == snpeffcalls[3]:
       if armname not in modifiercalls:
         modifiercalls[armname] = []
       modifiercalls[armname].append(pos)
     else: # was something other than "modifier"
       if armname not in notmodifiercalls:
         notmodifiercalls[armname] = []
       notmodifiercalls[armname].append(pos)

   # determine distances
   r_dist,s_dist = calc_distances(muts,arm_boundaries)

   # accumulate these results into grand tallies
   real_distances = [sum(x) for x in zip(r_dist,real_distances)]
   sim_distances = [sum(x) for x in zip(s_dist,sim_distances)]

   # calculate proportion of mutations in bins 50 or less
   p50 = calcperc(r_dist,50)
   real_50count.append([pid,sid,p50])
   p50 = calcperc(s_dist,50)
   sim_50count.append([pid,sid,p50])

   # tally modifier and non-modifier separately
   r_dist,s_dist = calc_distances(modifiercalls,arm_boundaries)
   mod_real_distances = [sum(x) for x in zip(r_dist,mod_real_distances)]
   mod_sim_distances = [sum(x) for x in zip(s_dist,mod_sim_distances)]
   r_dist,s_dist = calc_distances(notmodifiercalls,arm_boundaries)
   nmod_real_distances = [sum(x) for x in zip(r_dist,nmod_real_distances)]
   nmod_sim_distances = [sum(x) for x in zip(s_dist,nmod_sim_distances)]


#print "Real distances"
#print real_distances
#print "Simulated distances"
#print [x/float(nreps) for x in sim_distances]

outfile = open("simresults","w")
outline = "Distance\tRealdata\tSimdata\n"
outfile.write(outline)
freps = float(nreps)
for i in range(1,categories):
  outline = str(i) + "\t" + str(real_distances[i]) + "\t" + str(sim_distances[i]/freps) + "\n"
  outfile.write(outline)
outfile.close()

outfile = open("simresults_modifier","w")
outline = "ModDistance\tModRealdata\tModSimdata\n"
outfile.write(outline)
freps = float(nreps)
for i in range(1,categories):
  outline = str(i) + "\t" + str(mod_real_distances[i]) + "\t" + str(mod_sim_distances[i]/freps) + "\n"
  outfile.write(outline)
outfile.close()

outfile = open("simresults_notmodifier","w")
outline = "NotModDistance\tNotModRealdata\tNotModSimdata\n"
outfile.write(outline)
freps = float(nreps)
for i in range(1,categories):
  outline = str(i) + "\t" + str(nmod_real_distances[i]) + "\t" + str(nmod_sim_distances[i]/freps) + "\n"
  outfile.write(outline)
outfile.close()

outfile = open("simresults_50bin","w")
outline = "Pid\tSid\tReal\tSimulated\n"
outfile.write(outline)
for realval,simval in zip(real_50count,sim_50count):
  if realval[0] != simval[0]:
    print "disagreement on pid",realval[0],simval[0]
    outfile.close()
    exit()
  if realval[1] != simval[1]:
    print "disagreement on sid",realval[1],simval[1]
    outfile.close()
    exit()
  outline = realval[0] + "\t" + realval[1] + "\t" + str(realval[2])
  outline += "\t" + str(simval[2]) + "\n"
  outfile.write(outline)
outfile.close()



####

Multiple phase solution:

phase 1:  real data only, store (pickle!):
   overall distances per arm per sample
   mod distances per arm per sample
   nmod distances per arm per sample

phase 2:
   for each of overall, mod, and nmod:
     unpickle the relevant phase 1
     for each sample:
       conduct nreps whole simulations
         score overall P50 and per-chromosome P50 
       compare real P50 and real per-chromosome P50 to simulations
       accumulate real and simulated counts into grand tallies
       report on this sample
     report on grand tallies
     del the phase 1 variable!
     
