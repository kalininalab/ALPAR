#!/usr/bin/env python3

from __future__ import print_function
import sys
import os
import pygtrie

fileext=[".fa",".fas",".fastq",".mfasta",".fasta",".fsa",".fna"]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def unify(split, taxa):
    ls=len(split)
    lt=len(taxa)
    if ls==(lt-ls):
        sorted_inv = sorted([item for item in taxa if item not in split])
        sorted_split=sorted(split)
        if sorted_inv<sorted_split:
            return sorted_inv
        else:
            return sorted_split
    if (lt-ls)<ls:
        inv = [item for item in taxa if item not in split]
        return sorted(inv)
    else:
        return sorted(split)

def readtaxa(filename):
    taxa=set()
    for line in (s.strip().split()[0] for s in open(filename)):
        # remove file extensions
       if os.path.splitext(line)[1] in [".gz",".gzip",".zip"]:
            line=os.path.basename(os.path.splitext(line)[0])
       if os.path.splitext(line)[1] in fileext:
            line=os.path.splitext(line)[0]
       taxa.add(line)
    return(taxa)


def readfile(filename,taxa,min_size):
    splits=pygtrie.StringTrie()
    printmin=True
    for line in (s.strip() for s in open(filename)):
        #split line
        fields = line.split('\t')
        weight = float(fields[0])
        split = fields[1:]
        #store tax names if necessary
        for idx,f in enumerate(split):
            # remove file extensions
            if os.path.splitext(f)[1] in [".gz",".gzip",".zip"]:
                f=os.path.basename(os.path.splitext(f)[0])
            if os.path.splitext(f)[1] in fileext:
                f=os.path.splitext(f)[0]
                split[idx]=f
            if f not in taxa:
                eprint("taxa from split not in taxa: "+f)
                exit(1)
        split=unify(split,taxa)
        if  len(split)<min_size:
            if len(split)==0:
                eprint("ignored split \"all vs. none\"")
                continue
            if printmin: #print only on first occurrence
                eprint("ignored split smaller than specified minimum size of ", min_size)
                printmin=False
            continue
        split_txt="/".join(split)
        #update or store split
        if split_txt in splits:
            splits[split_txt]+=weight
        else:
            splits[split_txt]=weight

    return splits




def get_top(splitlist,l,num): #l=len(splitlist)
    

    # too short to be filtered at all?
    if l<=num:
        return splitlist
    # else..
    
    #sort by weights
    both_sorted=sorted([(w,s) for (s,w) in splitlist.items()], reverse=True)
    
    #cut top num
    for (_,s) in both_sorted[num:]:
        del splitlist[s]



def remove_leaves(splitlist,taxa):
    for t in taxa:
        if t in splitlist:
            del splitlist[t]


def split_comp(donor,reference,weighted):
    
    #sort by weights
    donor_sorted=sorted([(w,s) for (s,w) in donor.items()], reverse=True)
    max_donor=donor_sorted[0][0]

    ref_list=reference.keys()
    n=len(ref_list)
    w_all_ref=sum([w for (_,w) in reference.items()])
    max_ref=max([w for (_,w) in reference.items()])
    num_all_ref=len(reference.items())

    w_all=0.
    w_corr=0.
    w_corr_ref=0.
    num_all=0
    num_corr=0
    w_precision=-1.
    w_recall=-1.
    dist=0.
    w_dist=0.
    branch_score=0.


    if not weighted:
        print("#precision\trecall\tF1-score\tsymmetric_distance\tRF-distance\t(unweighted)")
    else:
        print("#precision\trecall\tF1-score\tsymmetric_distance\tbranch score\t(weighted)")

    for (w,s) in donor_sorted:
        w_all+=w
        num_all+=1
        w_ref=0.
        if s in ref_list:
            w_ref=reference[s]
            ref_list.remove(s)
            w_corr+=w
            num_corr+=1
            w_corr_ref+=reference[s]
 
        precision=(1.0*num_corr)/(1.0*num_all)
        recall=(1.0*num_corr)/(1.0*n)
        dist=((1.0*num_all-num_corr)/(1.0*num_all)) + ((1.0*num_all_ref-num_corr)/(1.0*num_all_ref))
        branch_score+=abs(w_ref/max_ref-w/max_donor)

        if w_all > 0:
            w_precision=(1.0*w_corr)/(1.0*w_all)
        if w_all_ref > 0 :
            w_recall=(1.0*w_corr_ref)/(1.0*w_all_ref)
            if w_all > 0:
                w_dist=((w_all-w_corr)/(w_all)) + ((w_all_ref-w_corr_ref)/(w_all_ref))

        if not weighted:
            print("\t".join([str(precision),str(recall),"-" if (precision*recall==0) else str(2*precision*recall/(precision+recall)),str(dist),str((num_all-num_corr)+(num_all_ref-num_corr))]))
        else: # weighted
            print("\t".join([str(w_precision),str(w_recall),"-" if (precision*recall==0) else str(2*w_precision*w_recall/(w_precision+w_recall)),str(w_dist),str(branch_score+sum([reference[s] for s in ref_list])/max_ref)]))


    eprint("\t".join([str(precision),str(recall),"-" if (precision*recall==0) else str(2*precision*recall/(precision+recall)),str(dist),str((num_all-num_corr)+(num_all_ref-num_corr)),"unweighted"]))
    eprint("\t".join([str(w_precision),str(w_recall),"-" if (precision*recall==0) else str(2*w_precision*w_recall/(w_precision+w_recall)),str(w_dist),str(branch_score+sum([reference[s] for s in ref_list])/max_ref),"weighted"]))
                     

                        

   
weighted=False
min_size=1


if len(sys.argv)>4:
	if sys.argv[4]=="-w":
		weighted=True
	else:
		min_size=int(sys.argv[4])
	if len(sys.argv)>5:
		if sys.argv[5]=="-w":
			weighted=True
		else:
			min_size=int(sys.argv[5])

if len(sys.argv)<4:
    eprint("Usage: comp.py <list of splits 1> <list of splits 2> <list of genomes (file names)>  [minimum split size] -w")
    eprint("list of splits in SANS output format: Per split one line: weight genomeA genomeB ... (fasta file extensions are ignored: "+", ".join(fileext)+")")
    eprint("Output on stdout: precision and recall of splits 1 w.r.t. splits 2 in terms of topological RF-distance; symmetric distance; Robinson-Foulds-Distance (unweighted) or branch score (weighted); for increasing number of splits considered, i.e., precision, recall and distance for *all* splits are in last line.")
    eprint("Only splits of size at least <minimum split size> are considered (default = 1, i.e. all; choose 2 to ignore trivial splits, i.e. leaf edges).")
    eprint("Default:     Precision = (number of splits 1 that are also in splits 2) / (total number of splits 1)")
    eprint("             Recall    = (number of splits 2 that are also in splits 1) / (total number of splits 2)")
    eprint("             F1-score  = Harmonic mean of precision and recall")
    eprint("             Distance  = (number of splits 1 that are not in splits 2) / (total number of splits 1)")
    eprint("                       + (number of splits 2 that are not in splits 1) / (total number of splits 2)") 
    eprint("             RF-Dist.  = number of splits 1 that are not in splits 2")
    eprint("                       + number of splits 2 that are not in splits 1")
    eprint("-w Weighted: Precision = (total weight of splits 1 that are also in splits 2) / (total weight of all splits 1)")
    eprint("             Recall    = (total weight of splits 2 that are also in splits 1) / (total weight of all splits 2)")
    eprint("             F1-score  = Harmonic mean of weighted precision and weighted recall")
    eprint("             Distance  = (total weight of splits 1 that are not in splits 2) / (total weight of splits 1)")
    eprint("                       + (total weight of splits 2 that are not in splits 1) / (total weight of splits 2)")
    eprint("             Branch score = Sum of weight differences per split")
    eprint("                            (Weights are normalized by maximum weight per input file, splits not occurring have weight zero.")
    eprint("Both unweighed and weighted statistics are output on stderr in any case.")
    sys.exit(1)


   
eprint("read taxa")
taxa=readtaxa(sys.argv[3])
eprint(len(taxa))

eprint("read split file 1")
splitlist1=readfile(sys.argv[1],taxa,min_size)
eprint("read split file 2")
splitlist2=readfile(sys.argv[2],taxa,min_size)
eprint("found:")
l1=len(splitlist1)
l2=len(splitlist2)
eprint(l1,l2)

#eprint("get copy without leaves")
#leaflesssplits1=splitlist1.copy()
#remove_leaves(leaflesssplits1,taxa)
#leaflesssplits2=splitlist2.copy()
#remove_leaves(leaflesssplits2,taxa)
#eprint("left:")
#ll1=len(leaflesssplits1)
#ll2=len(leaflesssplits2)
#eprint(ll1,ll2)

#eprint("filter top n**2-n ("+str(len(taxa)**2-len(taxa))+") from leafless sets")
#get_top(leaflesssplits1,ll1,len(taxa)**2-len(taxa))
#get_top(leaflesssplits2,ll2,len(taxa)**2-len(taxa))

#eprint("compute RF on leafless sets")
#rf_l=comp_diff_sqr(leaflesssplits1,leaflesssplits2,splits2lists(leaflesssplits1,leaflesssplits2))
#eprint("compute FM on leafless sets")
#fm_l=comp_diff_sqr(leaflesssplits1,leaflesssplits2,get_pairwise_ds(leaflesssplits1,leaflesssplits2,taxa))


#eprint("filter top n**2 ("+str(len(taxa)**2)+") from complete sets")
#get_top(splitlist1,l1,len(taxa)**2)
#get_top(splitlist2,l2,len(taxa)**2)

#eprint("compute RF")
#rf=comp_diff_sqr(splitlist1,splitlist2,splits2lists(splitlist1,splitlist2))
#eprint("compute FM")
#fm=comp_diff_sqr(splitlist1,splitlist2,get_pairwise_ds(splitlist1,splitlist2,taxa))
eprint("compute precision and recall:")
split_comp(splitlist1,splitlist2,weighted)


# output
#print("\t".join(map(str,[rf,rf_l,fm,fm_l,sc])))
#print("\t".join([str(diff_taxa)]))


