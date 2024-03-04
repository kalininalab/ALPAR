#!/usr/bin/env python3

from __future__ import print_function
import sys
import os


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if len(sys.argv)<3:
    eprint("Usage: sans2pdf.py <sans output file> <list of genomes (file names)>")
    eprint("Output: ")
    eprint("  NEXUS file <sans output file>.nexus")
    eprint("  PDF file <sans output file>.nexus.pdf")
    eprint("Requires:")
    eprint("  SplitsTree in your PATH")
    eprint("  sans2nexus.py in this directory")
    sys.exit(1)

nex=sys.argv[1]+".nexus"
pdf=nex+".pdf"

dir_path = os.path.dirname(os.path.realpath(__file__))

# sans -> nexus
os.system(" ".join([dir_path+"/sans2nexus.py",sys.argv[1],sys.argv[2],">",nex]))

# prepare SplitsTree command
f = open("splitstreecommands.tmp", "w")
print("OPEN FILE="+nex,file=f)
print("UPDATE",file=f)
print("UPDATE",file=f)
print("UPDATE",file=f)
print("EXPORTGRAPHICS format=PDF TEXTASSHAPES=YES file="+pdf+" REPLACE=yes",file=f)
print("QUIT",file=f)
f.close()

# nexus -> pdf
os.system("SplitsTree -g -c splitstreecommands.tmp")

# tidy up
os.system("rm -f splitstreecommands.tmp")
