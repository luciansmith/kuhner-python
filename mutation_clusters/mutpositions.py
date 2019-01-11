# mutposition.py          Mary Kuhner 2018/04/10

# choose 1 random sample from each patient, excluding gastric
# tabulate all mutation positions (2+ callers) for all those samples
# do NOT cast out duplicates
# pickle by chromosome


#####################################################################33
# main program

# where to find inputs
boundfile = "chrom_boundaries"
mutlocation = "/home/mkkuhner/seagate/data/wgs_dna/union/"

import random
import os

allowed_chroms = [str(x) for x in range(1,23)]
snpeffcalls = ["HIGH","MODERATE","LOW","MODIFIER"]


# read chromosome boundaries, needed to classify mutations by arm
boundaries = {}
arm_boundaries = {}
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


# locate all mutation files
allfiles = []
for root, dirs, files in os.walk(mutlocation):
  for file in files:
    if file.endswith(".annotated.txt"):
      allfiles.append(file)


# sort by patients and choose one file per patient at random
files_patient = {}
for file in allfiles:
  pid = file.split("-")[0]
  if pid not in files_patient:
    files_patient[pid] = []
  files_patient[pid].append(file)

chosen_files = []
for pid in files_patient:
  chosen_files.append(random.choice(files_patient[pid]))


# set up storage for results
mutstore_overall = {}
mutstore_mod = {}
mutstore_nmod = {}

for arm in arm_boundaries:
  mutstore_overall[arm] = []
  mutstore_mod[arm] = []
  mutstore_nmod[arm] = []


# read mutations from files
print "Assessing mutations in",len(chosen_files),"input files"
for mutfile in chosen_files:   
   muts = {}
   prefix = mutfile.split("/")[-1]
   prefix = prefix.split("-")
   sid = prefix[1]

   if sid.endswith("N"):
     continue
  
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

     mutstore_overall[armname].append(pos)

     snpeffcall = line[snpeffindex]
     if snpeffcall not in snpeffcalls:
       print "\nfound an abberant snpeffcall for chomosome",chr,
       print "position",pos,":  ",snpeffcall
       exit()
     if snpeffcall == snpeffcalls[3]:     # was MODIFIER
       mutstore_mod[armname].append(pos)
     else:                                # was not MODIFIER
       mutstore_nmod[armname].append(pos)

total_overall = 0
total_mod = 0
total_nmod = 0
for arm in mutstore_overall:
  total_overall += len(mutstore_overall[arm])
  total_mod += len(mutstore_mod[arm])
  total_nmod += len(mutstore_nmod[arm])
print "Found",total_overall,"mutations:",total_mod,"MOD and",total_nmod,"not MOD"
  
import pickle
pfile = "mutpositions.pkl"
picklefile = open(pfile, "w")
pickle.dump(mutstore_overall,picklefile)
pickle.dump(mutstore_mod,picklefile)
pickle.dump(mutstore_nmod,picklefile)
picklefile.close()
