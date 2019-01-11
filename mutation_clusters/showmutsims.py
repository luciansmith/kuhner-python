import pickle
nreps = 1000

pfile = open("indexscores.pkl","r")
overall_scores = pickle.load(pfile)
arm_scores = pickle.load(pfile)
overall_scores_mod = pickle.load(pfile)
arm_scores_mod = pickle.load(pfile)
overall_scores_nmod = pickle.load(pfile)
arm_scores_nmod = pickle.load(pfile)
pfile.close()

numsamples = len(overall_scores)

cutoff = 0.999
cutoff2 = 0.95
opval = 0
opval2 = 0
mpval = 0
mpval2 = 0
npval = 0
npval2 = 0
for score in overall_scores:
  if float(score)/nreps > cutoff:
    opval += 1
  if float(score)/nreps > cutoff2:
    opval2 += 1
for score in overall_scores_mod:
  if float(score)/nreps > cutoff:
    mpval += 1
  if float(score)/nreps > cutoff2:
    mpval2 += 1
for score in overall_scores_nmod:
  if float(score)/nreps > cutoff:
    npval += 1
  if float(score)/nreps > cutoff2:
    npval2 += 1

print "Results for p<",1-cutoff
print "Overall:",opval
print "Modifier:",mpval
print "Non-modifier:",npval

print "Results for p<",1-cutoff2
print "Overall:",opval2
print "Modifier:",mpval2
print "Non-modifier:",npval2

print "Arm results:  arm, p<",1-cutoff,", p<",1-cutoff2
for arm in arm_scores:
  opval = 0
  opval2 = 0
  for score in arm_scores[arm]:
    if float(score)/nreps > cutoff:
      opval += 1
    if float(score)/nreps > cutoff2:
      opval2 += 1
  print arm, float(opval)/numsamples, float(opval2)/numsamples

#import matplotlib.pyplot as plt
#
#figno = 0
#plt.figure(figno)
#plt.title("Overall results")
#plt.hist(overall_scores,bins=30,range=[0,nreps])
#
#figno = 1
#plt.figure(figno)
#plt.title("Overall results--Modifier")
#plt.hist(overall_scores_mod,range=[0,nreps],bins=30)
#
#figno = 2
#plt.figure(figno)
#plt.title("Overall results--Non-Modifier")
#plt.hist(overall_scores_nmod,range=[0,nreps],bins=30)
#
#
##for arm in arm_scores:
##  figno += 1
##  plt.figure(figno)
##  plt.title("Chrom"+arm)
##  plt.hist(arm_scores[arm])
#  
#plt.show()
