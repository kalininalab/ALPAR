#!/usr/bin/env python3


from __future__ import print_function
import sys
import os
import re
from math import log
from math import sqrt

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



def readfile(filename):
    taxa=False
    splits=False
    splitlist=list()
    idlist=dict()
    p_tax = re.compile('\[(\d+)\]\s\'(.+)\'')
    #[1, size=4]	65447.0	4 8 9 12,
    p_split= re.compile('\[.+\]\s+(\S+)\s+(.*),')
    for line in (s.strip() for s in open(filename)):
        #skip until "TAXLABELS"
        if not taxa and line=="TAXLABELS":
            taxa=True;
            continue;
        elif taxa:
            m = p_tax.search(line)
            if not m:
                taxa=False
                eprint("Taxa read:")
                eprint(len(idlist))
                continue
            else:
                idlist[int(m.group(1))]=m.group(2)
        elif not splits and line=="MATRIX":
            splits=True
            continue;
        elif splits:
            m= p_split.search(line)
            if not m:
                eprint("Splits read:")                
                eprint(len(splitlist))
                return splitlist
            else:
                splitlist.append(([idlist[int(i)] for i in m.group(2).split(' ')],m.group(1)))
    return None


if len(sys.argv)<2:
    eprint("Usage: nexus2sans.py <NEXUS file>")
    eprint("Output in SANS output format on stdout.")
    sys.exit(1)


eprint("read nexus file")
splitlist=readfile(sys.argv[1])


# output as tab sep file
for (s,w) in splitlist:
    print(w+"\t"+"\t".join(map(str,s)))

