#!/usr/bin/env python3

import dendropy
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_clusters(clusters):
	n=0
	for c in clusters:
		#if len(set([i.taxon.__str__().split('_')[1] for i in c]))>1:
			#eprint(set([i.taxon.__str__().split('_')[1] for i in c]))
		for t in [i.taxon.__str__()[1:len(i.taxon.__str__())-1] for i in c]:
			print(t+"\t"+str(n))
		n+=1

def eprint_cluster_hist(clusters,gold_standard):
	hist=dict()
	for c in clusters:
		#l=len(set([i.taxon.__str__().split('_')[1] for i in c]))
		l=len(set([gold_standard[i.taxon.__str__()[1:len(i.taxon.__str__())-1]] for i in c]))
		if not l in hist:
			hist[l]=1
		else:
			hist[l]+=1
		#if l==2:
			#for t in set([i.taxon.__str__().split('_')[1][0:len(i.taxon.__str__().split('_')[1])-1] for i in c]):
				#print(t)

	eprint("N_CL\tN_GOLDSTANDARD_CL")
	for k in sorted(hist.keys()):
		eprint(str(hist[k])+"\t"+str(k))

def eprint_bin_hist(clusters,gold_standard):
	clusters_per_bin=dict()
	for c in clusters:
		#for bin in set([i.taxon.__str__().split('_')[1] for i in c]):
		for bin in set([gold_standard[i.taxon.__str__()[1:len(i.taxon.__str__())-1]] for i in c]):
			if not bin in clusters_per_bin:
				clusters_per_bin[bin]=1
			else:
				clusters_per_bin[bin]+=1
	hist=dict()
	for bin in clusters_per_bin.keys():
		n=clusters_per_bin[bin]
		if not n in hist:
			hist[n]=1
		else:
			hist[n]+=1
	tot=0
	eprint("N_GOLDSTANDARD_CL\tN_CL")
	for k in sorted(hist.keys()):
		eprint(str(hist[k])+"\t"+str(k))
		tot+=hist[k]
	eprint("total number of gold standard clusters: "+ str(tot))


def proc_node(node,cluster,prev_edge_length):
	subclusters=[]
	# special case: just 0-edges and one other -> ignore and consider merged edge length
	#
	# ---------------
	#   a  \0   b
	#
	# =>
	#
	# ---------------
	#     a + b

	# count none-zero edges
	n=0
	nz_child=None
	for child in node.child_node_iter():
		if child.edge_length>0:
			nz_child=child
			n+=1

	# if only one, special treatment
	if n==1:
		# consider total weight
		nz_child.edge_length+=node.edge_length
		# continue with current cluster
		return proc_node(nz_child,cluster,prev_edge_length) 


	# which clusters do we get from below?
	for child in node.child_node_iter():
		subclusters.extend(proc_node(child,set(child.leaf_nodes()),node.edge_length))
	
	# is the edge leading to this node a split?
	if node.edge_length>=prev_edge_length:
		# yes, split:
		#  new subcluster
		#  which clusters do we get from below?
		#  remove those taxa from this subcluster
		for sc in subclusters:
			cluster.difference_update(sc)
		#  add remaining taxa as subcluster
		if len(cluster)>0:
			subclusters.append(cluster)
	
	return subclusters

def proc_root_roughly(root,cluster):
	subclusters=[]
	#chop each subtree of root as cluster
	for child in root.child_node_iter():
		subclusters.append(set(child.leaf_nodes()))
	return subclusters


# USAGE
if len(sys.argv)<2 or (len(sys.argv)>3 and "-r" not in sys.argv):
        eprint("Usage: newick2clusters.py <newick file> [-r] [<gold standard>]")
        eprint(" -r: more rough clustering")
        eprint(" gold standard: tab-separated file: taxon_name<tab>cluster_ID")
        eprint("                (output statistics on stderr if given)")
        eprint(" output: tab separated on stdout: taxon_name<tab>cluster_ID")
        sys.exit(1)


# read tree file
tree = dendropy.Tree.get(
        path=sys.argv[1],
        schema="newick",
        preserve_underscores=True)#,
        #suppress_internal_node_taxa=True,
        #suppress_leaf_node_taxa=True)

eprint("number of leaves: "+ str(len(tree.leaf_nodes())))



#set all none-length edges to 0
for edge in tree.postorder_edge_iter():
	if edge.length is None:
		edge.length = 0

#reroot to max degree node
root=tree.seed_node
for node in tree.preorder_node_iter():
	if node.num_child_nodes()>root.num_child_nodes():
		root=node
tree.reroot_at_node(root, update_bipartitions=True)

# CLUSTERING

if "-r" not in sys.argv:
	#compute clusters recursively
        clusters=proc_node(root,set(root.leaf_nodes()),0)
else: #rough
	#split below root
	clusters=proc_root_roughly(root,set(root.leaf_nodes()))
	
# OUTPUT
print_clusters(clusters)

# statistics
if (len(sys.argv)==3 and sys.argv[2]!="-r") or len(sys.argv)>3:
	#read file
	gold_standard=dict()
	if sys.argv[2]!="-r":
		fn=sys.argv[2]
	else:
		fn=sys.argv[3]
	for line in (s.strip() for s in open(fn)):
		#split line
		fields = line.split('\t')
		taxon = fields[0]
		cluster = fields[1]
		if taxon in gold_standard.keys():
			eprint("ERROR in reading gold standard. Taxon "+taxon+" is assigned twice.")
			exit(1)
		#elif not taxon in tree.leaf_nodes():
			#eprint("ERROR in reading gold standard. Taxon "+taxon+" from gold standard not in newick file.")
			#exit(1)
		else:
			gold_standard[taxon]=cluster
	#check if all taxons defined
	#eprint(gold_standard.keys())
	for n in tree.leaf_nodes():
		t=n.taxon.__str__()[1:len(n.taxon.__str__())-1]
		if t not in gold_standard.keys():
			eprint("ERROR in reading gold standard. Taxon "+t+" from newick file not in gold standard file.")
			exit(1)

	eprint_cluster_hist(clusters,gold_standard)
	eprint_bin_hist(clusters,gold_standard)
