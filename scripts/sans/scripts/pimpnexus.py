#!/usr/bin/env python3


from __future__ import print_function
import sys
import os
import colorsys

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)


def readcolorfile(filename):
	label_map=dict()
	labels=[]
	for line in (s.strip() for s in open(filename)):
		#split line
		fields = line.split('\t')
		if len(fields)<2:
			fields.append("N/A")
		if not fields[1] in labels:
			labels.append(fields[1])
		label_map[fields[0]]=labels.index(fields[1])
	# map color id to rgb
	color_map=dict()
	N = len(labels)
	HSV_tuples = [(x/N, 1, 1) for x in range(N)]
	RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
	for key in label_map.keys():
		rgb=RGB_tuples[label_map[key]]
		color_map[key]=' '.join(list(map(lambda c: str(round(c*255)), rgb)))
		label_map[key]=labels[label_map[key]]
	return (color_map,label_map,labels)


#if len(sys.argv)>3:
#    eprint(">3")
#else:
#    eprint("<=3")





if len(sys.argv)<3:
    eprint("Usage: pimpnexus.py <nexus file> <leaf colors file>")
    eprint("leaf colors file: tab seperated: taxon, color (any kind of ID)")
    eprint("Output in NEXUS format on stdout.")
    eprint("If vertex is already colored, coloring won't be changed.")
    eprint("If several taxa of different color are at one leaf vertex, an empty circle is used as label.")
    sys.exit(1)


eprint("read color file")
(color_map, label_map, labels)=readcolorfile(sys.argv[2])
eprint(len(color_map)," color mappings read.")



eprint("pimp nexus file")
translate=False
vertices=False
vlabels=False
edges=False
taxa_map=dict()
minx=0
maxx=0
miny=0
maxy=0
maxvertexid=0
nvertices=0
for line in (s.strip() for s in open(sys.argv[1])):
	if line==";":
		#output legend
		if vertices:
			N = len(labels)
			HSV_tuples = [(x/N, 1, 1) for x in range(N)]
			RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
			for i in range(len(labels)):
				maxvertexid+=1
				rgb=RGB_tuples[i]
				print(str(maxvertexid)+" "+str(maxx)+" "+str(maxy+(i*(maxy-miny)/20))+" w=10 h=10 bg="+' '.join(list(map(lambda c: str(round(c*255)), rgb)))+" fg="+' '.join(list(map(lambda c: str(round(c*255)), rgb)))+",")
		if vlabels:
			for i in range(len(labels)):
				 #l=9 x=26 y=-25 f='Dialog-PLAIN-10',
				print(str(nvertices+i+1)+" '"+labels[i]+"' l=9 x=15 y=5 f='Dialog-PLAIN-10',")
		# end legend
		translate=False
		vertices=False
		vlabels=False
	
	if translate:
		fields = line[:len(line)-1].split(' ')
		vid=fields[0]
		taxa=list(map(lambda s: s[1:len(s)-1], fields[1:]))
		taxa_map[vid]=taxa
		print(line)
	elif vertices:
		maxvertexid+=1;
		fields = line[:len(line)-1].split(' ')
		if float(fields[1])<minx:
			minx=float(fields[1])
		if float(fields[1])>maxx:
			maxx=float(fields[1])
		if float(fields[2])<miny:
			miny=float(fields[2])
		if float(fields[2])>maxy:
			maxy=float(fields[2])
		if line.find("bg")==-1 and line.find("fg")==-1 and fields[0] in taxa_map.keys():
			if taxa_map[fields[0]][0] in color_map.keys() and len(set(map(lambda x: color_map[x], taxa_map[fields[0]])))==1:
				print(fields[0]+" "+fields[1]+" "+fields[2]+" w=10 h=10 bg="+color_map[taxa_map[fields[0]][0]]+" fg="+color_map[taxa_map[fields[0]][0]]+",")
			else:
				print(fields[0]+" "+fields[1]+" "+fields[2]+" w=10 h=10 fg=0 0 0 bg=255 255 255,")
		else:
			print(line)
	elif vlabels:
		dummy=0
		#fields = line[:len(line)-1].split(' ')
		#tax=fields[0]
		#if taxa_map[tax][0] in color_map.keys():
			#print(line[:len(line)-1]+" lc='"+color_map[taxa_map[tax][0]]+"',")
		#else:
			#print(line)
	elif edges:
		if line.find("w=")==-1:
			print(line)
		else:
			print(line[:line.find("w=")-1]+",")
			
	elif line.find("DIMENSIONS")>-1 and line.find("nvertices=")>-1:
		#DIMENSIONS ntax=239 nvertices=1033 nedges=1599;
		prefix=line.split("nvertices=")[0]
		nvertices=int(line.split("nvertices=")[1].split(" ")[0])
		suffix=line.split("nvertices=")[1].split(" ")[1]
		print(prefix+"nvertices="+str(nvertices+len(labels))+" "+suffix)
	else:
		print(line)
	if line=="TRANSLATE":
		translate=True
	if line=="VERTICES":
		vertices=True
	if line=="VLABELS":
		vlabels=True
	if line=="EDGES":
		edges=True
