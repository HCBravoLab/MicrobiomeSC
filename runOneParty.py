import sys
import subprocess

f = open(sys.argv[1])

line = f.readlines()


subprocess.call(line[(int)(sys.argv[2])].split("&")[0], shell=True)
