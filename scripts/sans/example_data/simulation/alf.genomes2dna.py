#!/usr/bin/env python3

from sys import stdout,stderr,exit,argv
from optparse import OptionParser
# from random import gammavariate as gamma, choice 
from os.path import basename
# from cStringIO import StringIO
import logging

from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import re

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

FASTA_HEADER_PAT = re.compile('^(G\d+_SE\d+), sequence type: (.*), locus: (-?\d+)$')

DEFAULT_MIN_LENGTH = 100
DEFAULT_CHROMOSOME_NO = 6

if __name__ == '__main__':

    usage = 'usage: %prog [options] <{ALF SIMULATED GENOME}_dna.fa>'
    parser = OptionParser(usage=usage)
#    parser.add_option('-i', '--insert_random_dna', dest='insertRandom',
#            default=-1, type=int, 
#            help='Insert random DNA between markers of minimum length <i> ' + \
#                    'as to simulate highly diverged non-coding areas. ' + \
#                    'Negative values for <i> deactivate this option. ' + \
#                    '[default: %default]')
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    
    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    
    cf = logging.FileHandler('%s.log' %(basename(argv[0]).rsplit('.', 1)[0]), mode='w')
    cf.setLevel(logging.INFO)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t++ %(message)s'))

    LOG.addHandler(cf)
    LOG.addHandler(ch)

    #
    # main 
    #

    
    genome = list()

    for rec in SeqIO.parse(open(args[0]), 'fasta'):
        m = FASTA_HEADER_PAT.match(rec.description)
        if not m:
            print >> stderr, 'Unable to parse FASTA sequence header "%s". Exiting' %(rec.description)
            exit(1)
        locus = int(m.group(3))
        if locus < 0:
            locus = -locus
            genome.append((locus, str(rec.seq.reverse_complement())))
        else:
            genome.append((locus, str(rec.seq)))

    genome.sort()
    _, dna = zip(*genome) 
#    if options.insertRandom >= 0:
#        import pdb; pdb.set_trace() 
#        inter_marker_len = map(lambda x: x >= options.insertRandom and x or 0,
#                [int(gamma(0.7, 5000)) for _ in xrange(len(dna))])
#        random_dna  = lambda x: ''.join(choice('ACGT') for _ in xrange(x))
#        i = 0
#        for l in inter_marker_len:
#            if l:
#                dna.insert(i, random_dna(l))
#                i += 2
#            else:
#                i += 1

    gName = basename(args[0]).rsplit('_dna.fa', 1)[0]
    gRec = SeqRecord(Seq(''.join(dna), generic_dna), 
                id='%s|chromosome|1' %gName, description='')
    SeqIO.write(gRec, stdout, 'fasta')

