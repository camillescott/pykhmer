#!/usr/bin/env python

import screed
import argparse
import sys

from dna_sequence import DNASequence as dna
from kmer_cms import KMerCountMinSketch as kcms

def normalize_by_median(cm, fn, cutoff, o_fn):

    n_kept = 0
    with open(o_fn, 'wb') as outfp:
        gen = (record for record in screed.open(fn))
        for record in cm.normalize_by_median(gen, cutoff):
            n_kept += 1
            outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
    
    print >>sys.stderr, 'Kept', n_kept, 'sequences'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_file')
    parser.add_argument('-o', dest='output_file')
    parser.add_argument('-C', dest='cutoff', default=10, type=int)
    parser.add_argument('-k', dest='ksize', default=20, type=int)
    parser.add_argument('-x', dest='htsize', default=1e6, type=float)
    parser.add_argument('-N', dest='nhashes', default=4, type=int)
    args = parser.parse_args()

    print 'Running normalize by median with', args.cutoff, args.ksize, args.htsize, args.nhashes, 'on', args.input_file

    cm = kcms(args.ksize, args.htsize, args.nhashes)

    if not args.output_file:
        args.output_file = args.input_file + '.keep'
    normalize_by_median(cm, args.input_file, args.cutoff, args.output_file)

if __name__ == '__main__':
    main()
