import sys
import subprocess

f = open(sys.argv[1])

line = f.readlines()

name = ["gen", "", "eva"]
p = (int)(sys.argv[2])
cmd = line[p].split("$")[0]+"| tee result/"+sys.argv[1].split(".")[0]+"_"+name[p]
print cmd
subprocess.call(cmd, shell=True)
