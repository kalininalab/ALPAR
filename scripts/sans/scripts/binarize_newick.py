#!/usr/bin/env python3

import dendropy
import sys






def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)




if len(sys.argv)<2:
    eprint("Usage: binarize_newick.py <newick file> ")
    eprint("Arbitrarily resolves polytomies using 0-length edges.")
    eprint("Output in newick format on stdout.")
    sys.exit(1)


tree = dendropy.Tree.get(
        path=sys.argv[1],
        schema="newick",
        preserve_underscores=True)

eprint("number of nodes read: "+ str(len(tree.leaf_nodes())))

 #resolve_polytomies(limit=2, update_bipartitions=False, rng=None)[source]

    #Arbitrarily resolve polytomies using 0-length edges.
    #Parameters:	

        #limit (int) – The maximum number of children a node can have before being resolved.
        #update_bipartitions (bool) – If True, then bipartitions will be calculated.
        #rng (random.Random object or None) – If rng is an object with a sample() method then the polytomy will be resolved by sequentially adding, generating all tree topologies equiprobably. rng.sample() should behave like random.sample() If rng is None, then polytomy is broken deterministically by repeatedly joining pairs of children.


tree.resolve_polytomies(2)

print(tree.as_string("newick"))
