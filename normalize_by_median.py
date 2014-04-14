#!/usr/bin/env python

import screed
import argparse
import sys

from kmer_cms import KMerCountMinSketch as kcms
from pykhmer_utils import iter_kmers

def normalize_by_median(cm, fn, cutoff):

    n_kept = 0
    total = 0
    with open(fn + '.keep', 'wb') as outfp:
        for n, record in enumerate(screed.open(fn)):
            if n % 10000 == 0:
                print >>sys.stderr, n, 'sequences processed'
            total += 1
            med = cm.get_median_count(record.sequence)
            if med > cutoff:
                continue
            else:
                n_kept += 1
                cm.count_sequence(record.sequence)
                outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
    print >>sys.stderr, 'Kept', n_kept, 'of', total, 'sequences'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_file')
    parser.add_argument('-C', dest='cutoff', default=10, type=int)
    parser.add_argument('-k', dest='ksize', default=20, type=int)
    parser.add_argument('-x', dest='htsize', default=1e6, type=float)
    parser.add_argument('-N', dest='nhashes', default=4, type=int)
    args = parser.parse_args()

    cm = kcms(args.ksize, args.htsize, args.nhashes)

    normalize_by_median(cm, args.input_file, args.cutoff)

if __name__ == '__main__':
    main()
