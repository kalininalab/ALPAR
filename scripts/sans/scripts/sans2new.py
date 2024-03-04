#!/usr/bin/env python3


from __future__ import print_function
import sys
import os
import pygtrie
from math import log
from math import sqrt

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if len(sys.argv)<2 or (len(sys.argv)>=3 and sys.argv[2]!="-g") or len(sys.argv)>4:
    eprint("Usage: sans2new.py <SANS file> [-g]")
    eprint("SANS file = tab separated list of splits: <weight> <tax1> <tax2> ...")
    eprint("-g [<filename>] greedy filter compatible subset of splits [and write splits to <filename>]")
    eprint("Output in NEWICK format on stdout.")
    sys.exit(1)



def readfile(filename):
    taxa=dict()
    inv_taxa=dict()
    i=1
    splits=pygtrie.StringTrie()
    for line in (s.strip() for s in open(filename)):
        #split line
        fields = line.split('\t')
        weight = float(fields[0])
        split_orig = fields[1:]
        #store tax names if necessary
        split=[] #storing ids
        for f in split_orig:
            if f not in taxa.keys():
                taxa[f]=i
                inv_taxa[i]=f
                split.append(i)
                i+=1
            else:
                split.append(taxa[f])
        split_txt="/".join(map(str,split))
        #update or store split
        if split_txt in splits:
            splits[split_txt]+=weight
        else:
            splits[split_txt]=weight
 
    return (splits,inv_taxa)

def findInTree(currentset,split,allTaxa):
    (subsets,taxa)=currentset
    # possible cases:
    # splitsize <2: nothing to be done
    # split equals one subset -> warning: split twice
    # split is fully contained in one subset -> recurse
    # inverse split ... (i.e. split covers one subset partially) -> recurse with inverse
    # split covers several subsets completely -> introduce new split
    
    if len(split)<2 or len(allTaxa)-len(split)<2:
        return True
    
    fullycoveredsubsets=[]
    partiallycoveredsubset=None
    if len(subsets)==0: # leaf set
            eprint("ERROR: leaf = split should be ignore in the first place!?")
            exit(1)
    for subset in subsets:
        (_,subtaxa)=subset
        if split == subtaxa:
            eprint("WARING: split observed twice: ",split)
            return True
        if split.issubset(subtaxa):
            return findInTree(subset,split,allTaxa)
        if subtaxa.issubset(split):
            fullycoveredsubsets.append(subset)
        elif not subtaxa.isdisjoint(split): # does intersect
            if partiallycoveredsubset:
                return False # there cannot be more than one
            else:
                partiallycoveredsubset=subset
                
    if partiallycoveredsubset:
        if len(fullycoveredsubsets)==len(subsets)-1:
            #recurse into this subset with inverse split
            inversesplit=allTaxa.difference(split)
            if inversesplit.issubset(partiallycoveredsubset[1]):
                return findInTree(partiallycoveredsubset,inversesplit,allTaxa)
            else:
                return False
        else:
            return False
    elif len(fullycoveredsubsets)>1:
        # introduce new split
        newsubset=fullycoveredsubsets
        newsubtaxa=set()
        for (_,subtaxa) in fullycoveredsubsets:
            newsubtaxa.update(subtaxa)
        newset=(newsubset,newsubtaxa)
        # remove old sets and add new set
        for subset in fullycoveredsubsets:
            subsets.remove(subset)
        subsets.append(newset)
        return True
    else:
        eprint("this cannot be: just one fully covered subset and nothing else!?")
        eprint(split)
        eprint(fullycoveredsubsets[0])
        eprint(subsets)
        exit(1)
        
    eprint("ERROR: This point in code should never be reached")
    exit(1)
    
def printtree(currentset,allTaxa):
    (subsets,taxa)=currentset
    if len(subsets)==0: # leaf set
        if len(taxa)==0:
            eprint("error: child with no taxon!?")
            exit(1)
        elif len(taxa)==1:
            if allTaxa:
                return allTaxa[list(taxa)[0]]
            else:
                return str(list(taxa)[0])
        else:
            eprint("error: child with more than one taxon!?")
            exit(1)
    else:
        s="("
        s+=",".join(map(lambda ss: printtree(ss,allTaxa), subsets))
        s+=")"
        return s

# Read input
eprint("read split file")
(splitlist,alltaxa)=readfile(sys.argv[1])
eprint("read:")
l=len(splitlist)
eprint(l, "splits")
eprint(len(alltaxa), "taxa")

# set: (list of subsets, set of all contained taxa)
# leaf set: (empty list, set on one taxon)
#complete set
subsets=[]
for taxon in alltaxa.keys():
    newset=([],set([taxon]))
    subsets.append(newset)
taxa=set(alltaxa.keys())
sets=(subsets,taxa)


if len(sys.argv)>=3 and sys.argv[2]=="-g": #greedy filter
    
    splitfile=None
    if len(sys.argv)>=4:
            splitfile=open(sys.argv[3],"w+")
    
    # sort given splits by weight
    #split initial set by splits
    
    
    weightedsplits= sorted([(w,s) for (s,w) in splitlist.items()], reverse=True)
    splitcount=0
    for (w,s) in weightedsplits:
        split=set(map(lambda x: int(x), s.split('/')))
        #split if possible
        if findInTree(sets,split,set(alltaxa.keys())):
            splitcount+=1
            if splitfile:
                print(str(w)+"\t"+"\t".join(map(lambda t: alltaxa[t],split)),file=splitfile)
    if splitfile:
            splitfile.close()
    eprint(splitcount, " compatible splits (of ",(2*len(alltaxa)-3)," possible splits in a fully resolved tree)")
    s=printtree(sets,alltaxa)
    s+=";"
    print(s)
    
else: # transform to newick and exit if incompatible
    
    #split initial set by splits
    splits=[set(map(lambda x: int(x), s)) for s in [s.split('/') for (s,_) in splitlist.items()]]

    for split in splits:
        #split if possible
        if not findInTree(sets,split,set(alltaxa.keys())):
            eprint("Splits are incompatible!")
            eprint("Tree so far: ",printtree(sets,alltaxa))
            eprint("Encountered incompatible split: ", set([alltaxa[t] for t in split]))
            exit(1)
            
    eprint("compatible")
    s=printtree(sets,alltaxa)
    s+=";"
    print(s)
    

