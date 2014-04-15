#!/usr/bin/env python

import array
import sys

from dna_sequence import DNASequence as dna
from pykhmer_utils import get_n_primes_above_x as get_primes
#from kmer_hash import hash_forward, hash_reverse, smhash_forward
from smhasher import murmur3_x86_64 as smhash

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

        kmer_hashed = smhash(kmer)

        return min(self.count_arrays[i][kmer_hashed % p] for i, p in enumerate(self.tablesizes))

    def count(self, kmer):
        kmer_hashed = smhash(kmer)
        try:
            for array_id, table_size in enumerate(self.tablesizes):
                self.count_arrays[array_id][kmer_hashed % table_size] += 1        
        except IndexError as e:
            print >>sys.stderr, 'error couning kmer', kmer,
            print >>sys,stderr, e

    def count_sequence(self, sequence):
        count = self.count
        for n, kmer in enumerate(sequence.kmers(self.ksize)):
            count(kmer)
        return n+1

    def get_median_count(self, sequence):
        get = self.get
        counts = sorted(get(kmer) \
                    for kmer in sequence.kmers(self.ksize))
        return counts[len(counts) // 2]

    def normalize_by_median(self, seq_gen, cutoff):
        median = self.get_median_count
        count_sequence = self.count_sequence
        for record in seq_gen:
            seq = dna(record.sequence)
            if median(seq) < cutoff:
                yield record
            else:
                count_sequence(seq)
