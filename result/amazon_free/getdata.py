import os


name = ["alphaDiversity", "chiSquare", "differentialAbundance", "oddsRatio"]
data = ["HMP", "PGP", "MSD"]


print "Running Time"
print "\tHMP Generator\tHMP Evaluator\tPGP Generator\tPGP Evaluator\tMSD Generator\tMSD Evaluator"
for n in name:
   print n,"\t",
   for d in data:
      for line in open(n+d+"_gen").readlines():
         if "running" in line:
            print line.split(":")[1].strip(),"\t",
      for line in open(n+d+"_eva").readlines():
         if "running" in line:
            print line.split(":")[1].strip(),"\t",
   print ""

print "\n","Circuit Size"
print "\tHMP\tPGP\tMSD"
for n in name:
   print n,"\t",
   for d in data:
      for line in open(n+d+"_eva").readlines():
         if "Gates" in line:
            print line.split(":")[1].strip(),"\t",
   print ""

print "\n", "Network Bandwidth"
print "\tHMP to Generator\tHMP to evaluator\tPGP to Generator\tPGP to Evaluator\tMSD to Generator\tMSD to Evaluator"
for n in name:
   print n,"\t",
   for d in data:
      for line in open(n+d+"_eva").readlines():
         if "Data" in line:
            print line.split(":")[1].split("MB")[0].strip(),"\t",
   print ""

