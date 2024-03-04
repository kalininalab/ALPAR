#!/usr/bin/env python3

from __future__ import print_function
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
import sys
import os
from random import randint

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)


def replace(c):
	r=randint(0,2)
	if c.upper()=='A':
		return "RWM"[r]
	elif c.upper()=='C':
		return "YSM"[r]
	elif c.upper()=='G':
		return "RSK"[r]
	elif c.upper()=='T':
		return "YWK"[r]
	else:
		eprint("Unexpected character: "+c)
		exit(1)

r=0.001 #ratio of replced symbols
n=False #do replace by R,W,M,Y,S,K (False) or N (True)

if len(sys.argv)<2 or sys.argv[1] in ["h","-h","--help", "-help"]:
	eprint("Usage: IUPACize.py <multi fasta file> [ratio] [-n]")
	eprint("Replaces a DNA base by N with prob. <ratio>.")
	eprint("-n replace by N; otherwise: character representing a second base")
	eprint("Default ratio = "+str(r))
	eprint("Output on stdout")
	exit(1)

if len(sys.argv)>2:
	if sys.argv[2]=="-n":
		n=True
		if len(sys.argv)>3:
			r=float(sys.argv[3])
	else:
		r=float(sys.argv[2])
		if len(sys.argv)>3:
			if sys.argv[3]=="-n":
				n=True
			else:
				eprint("invalid parameter: "+sys.argv[3])
				exit(1)


eprint("Using ratio of "+str(r))

records=SeqIO.parse(sys.argv[1], "fasta")

#read and manipulate
for record in records:
	length=len(record)
	eprint("Processing "+record.id+" of length "+str(length))
	#modify sequence
	#how many chacacters?
	num=round(length*r)
	eprint("modifying "+str(num)+" characters")
	positions=set()
	while len(positions)<num:
		positions.add(randint(0,length-1))
	s=str(record.seq)
	if n:
		for pos in positions:
			s= s[:pos-1] + 'N' + s[pos:]
	else:
		for pos in positions:
			s= s[:pos-1] + replace(s[pos]) + s[pos:]
	#set new sequence
	record.seq=Seq(s)
	#output
	SeqIO.write(record,sys.stdout, "fasta")
