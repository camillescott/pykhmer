#!/usr/bin/env python

import array
import sys

from pykhmer_utils import get_n_primes_above_x as get_primes
from pykhmer_utils import iter_kmers
from kmer_hash import hash_forward, hash_reverse

class KMerException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return 'KMerException: ' + self.msg

class KMerCountMinSketch(object):

    def __init__(self, ksize, minhashsize, nhashes, count_type='H'):
        self.ksize = ksize
        self.minhashsize = long(minhashsize)
        self.nhashes = nhashes

        self.count_type = count_type

        self.tablesizes = get_primes(self.nhashes, self.minhashsize)
        print self.tablesizes

        self.count_arrays = []
        for table_size in self.tablesizes:
            self.count_arrays.append(array.array(count_type, [0] * long(table_size)))

    
    def get(self, kmer):
        if len(kmer) != self.ksize:
            raise KMerException('Tried to get count \
                                with invalid k-mer size')
        kmer_hashed,_ = hash_forward(kmer, self.ksize)
        counts = []
        try:
            for array_id, table_size in enumerate(self.tablesizes):
                counts.append(self.count_arrays[array_id][kmer_hashed % table_size])
        except IndexError as e:
            print >>sys.stderr, 'error counting kmer', kmer,
            print >>sys,stderr, e    

        return min(counts)

    def count(self, kmer):
        kmer_hashed,_ = hash_forward(kmer, self.ksize)
        try:
            for array_id, table_size in enumerate(self.tablesizes):
                self.count_arrays[array_id][kmer_hashed % table_size] += 1        
        except IndexError as e:
            print >>sys.stderr, 'error couning kmer', kmer,
            print >>sys,stderr, e

    def count_sequence(self, sequence):
        for n, kmer in enumerate(iter_kmers(sequence, self.ksize)):
            self.count(kmer)
        return n

    def get_median_count(self, sequence):
        counts = [self.count(kmer) \
                    for kmer in iter_kmers(sequence, self.ksize)]
        counts.sort()
        return counts[len(counts) // 2]

    
    

