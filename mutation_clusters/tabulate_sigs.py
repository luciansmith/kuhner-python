# this program uses previous runs of deconstructsigs, plus info about
# the number of mutations provided by the deconstruct_input.py program,
# to make a report for Kenji in the following format:

# pid, sid, number clustered mutations, proportion in each signature

#################
# process numbers of mutations

numfile = "total_clustered.csv"

mutdict = {}
for line in open(numfile,"r"):
  line = line.rstrip().split(",")
  num = line[1]
  id = line[0].split("_")
  pid = id[0]
  sid = id[1]
  if sid.endswith("N"):  continue    # normal sample
  if pid not in mutdict:
    mutdict[pid] = {}
  mutdict[pid][sid] = num     # minus the header


##################
# process deconstructsigs results

outfilename = "clustersigs.csv"
outfile = open(outfilename,"w")

header = "PID_SID,nummuts"
for i in range(1,31):
  header += ",signature" + str(i)
header += "\n"
outfile.write(header)

pids = sorted(mutdict.keys())
for pid in pids:
  sids = sorted(mutdict[pid].keys())
  for sid in sids:
    outline = pid + "_" + sid + "," + mutdict[pid][sid]
    infilename = "signature_runs/" + pid + "_" + sid + "_sigs.tsv"
    signatures = open(infilename,"r").readlines()
    for i in [1, 3, 5, 7, 9, 11]:
      line = signatures[i].rstrip().split()
      for j in [1, 2, 3, 4, 5]:
        outline += "," + line[j] 
    outline += "\n"
    outfile.write(outline)

outfile.close()
      
