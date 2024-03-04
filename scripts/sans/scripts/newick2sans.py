#!/usr/bin/env python3

import dendropy
import os
import sys


fileext=[".fa",".fas",".fastq",".mfasta",".fasta",".fsa",".fna"]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



if len(sys.argv)<2:
    eprint("Usage: newick2sans.py <newick file> [<list of genomes (file names) to be considered, i.e., prune tree to given list>]")
    eprint("Output in SANS output format on stdout.")
    sys.exit(1)


tree = dendropy.Tree.get(
        path=sys.argv[1],
        schema="newick",
        preserve_underscores=True)

eprint("number of nodes read: "+ str(len(tree.leaf_nodes())))


if len(sys.argv)>2:
    eprint("pruning tree...")
    leaflables=[nd.taxon.label for nd in tree.leaf_nodes()]
    #read taxa
    taxa=set()
    for line in (s.strip() for s in open(sys.argv[2])):
        # remove file extensions
        if os.path.splitext(line)[1] in [".gz",".gzip",".zip"]:
            line=os.path.basename(os.path.splitext(line)[0])
        if os.path.splitext(line)[1] in fileext:
            line=os.path.splitext(line)[0]
        taxa.add(line)
        if line not in leaflables:
            eprint("taxa from list not in tree: "+line)
    eprint("taxa read: "+str(len(taxa)))
    node_filter_fn = lambda nd: nd.is_internal() or nd.taxon.label in taxa
    tree = tree.extract_tree(node_filter_fn=node_filter_fn)


tree.encode_bipartitions(suppress_unifurcations=True,collapse_unrooted_basal_bifurcation=True)
for edge in tree.preorder_edge_iter():
    # omit root edge -> split "all vs. none"
    if len(edge.bipartition.leafset_taxa(tree.taxon_namespace)) == len(tree.leaf_nodes()):
        continue;
    if edge.length:
        split=str(edge.length)+"\t"
    else:
        split="0.0\t"
    split+="\t".join([t.label for t in edge.bipartition.leafset_taxa(tree.taxon_namespace)])
    print(split)
