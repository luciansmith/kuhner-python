import os

# get list of Strelka files
targets = []
for root,dirs,files in os.walk("/home/mkkuhner/seagate/data/wgs_dna/filtered_strelka"):
  for file in files:
    if file.endswith("vcf"):
      targets.append(file)

# get list of BAM files
bamnames = []
for line in open("curated_bamlist"):
  bamnames.append(line.rstrip())

# pull BAM files corresponding to Strelka files
# watch out for 391-23521-155R-31

buckettargets = open("buckettargets.txt","w")
for target in targets:
  target = target.rstrip().split("-")
  # special case code for the pilot sample included in study
  if (target[0] == "391" and target[1] == "23521"):
    bam = "readonly/P01CA91955-WGS80/Project_REI_12051_11321_391_611_B01_SOM_WGS.2017-08-04/Sample_391-23521-155R-31/analysis/391-23521-155R-31.final.bam\n"
    buckettargets.write(bam)
    continue
  fullname = "-".join(target[0:4])
  found = False
  for bam in bamnames:
    if fullname in bam:
#      print("Found a match for",fullname)
#      print(bam)
      bam = bam + "\n"
      buckettargets.write(bam)
      found = True
      break
  if not found:
    print("************ No match for",fullname)

buckettargets.close()
