#!/usr/bin/env python

import array
import sys

from dna_sequence import DNASequence as dna
from pykhmer_utils import get_n_primes_above_x as get_primes
#from kmer_hash import hash_forward, hash_reverse, smhash_forward
from smhasher import murmur3_x64_128 as smhash

class KMerException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return 'KMerException: ' + self.msg

class KMerCountMinSketch(object):

    def __init__(self, ksize, minhashsize, nhashes, count_type='B'):
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

    def n_occupied(self):
        return min(sum(arr) for arr in self.count_arrays)

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

import numpy as np

class KMerCountMinSketch_numpy(object):

    def __init__(self, ksize, minhashsize, nhashes, dtype=np.uint8):
        self.ksize = ksize
        self.minhashsize = long(minhashsize)
        self.nhashes = nhashes

        self.dtype = dtype
    
        self.tablesizes = np.array(get_primes(self.nhashes, self.minhashsize), dtype=np.uint64)
        print self.tablesizes

        self.row_indices = np.arange(self.nhashes)
        self.counts = np.zeros((nhashes, max(self.tablesizes)), dtype=dtype)
        
    def get(self, kmer):
        if len(kmer) != self.ksize:
            raise KMerException('Tried to get count \
                                with invalid k-mer size')

        kmer_hashed = smhash(kmer)
        return self.counts[self.row_indices, np.mod(kmer_hashed, self.tablesizes)].min()

    def get_local(self, kmer):
        kmer_hashed = smhash(kmer)
        return min(self.counts[i, kmer_hashed % p] for i, p in enumerate(self.tablesizes))
     

    def count(self, kmer):
        kmer_hashed = smhash(kmer)
        kmer_bins = np.mod(kmer_hashed, self.tablesizes)
        self.counts[self.row_indices,kmer_bins] += 1

    def count_local(self, kmer):
        kmer_hashed = smhash(kmer)
        kmer_hashed = smhash(kmer)
        for i, p in enumerate(self.tablesizes):
            self.counts[i,kmer_hashed % p] += 1

    def count_sequence(self, sequence):
        count = self.count
        for n, kmer in enumerate(sequence.kmers(self.ksize)):
            count(kmer)
        return n+1

    def count_sequence_local(self, sequence):
        rows = self.row_indices
        ts = self.tablesizes
        for kmer in sequence.kmers(self.ksize):
            kh = smhash(kmer)
            self.counts[rows, kh % ts] += 1

    def n_occupied(self):
        return min(np.count_nonzero(row) for row in self.counts)

    def get_median_count(self, sequence):
        get = self.get
        return np.median(get(kmer) for kmer in sequence.kmers(self.ksize))

    def normalize_by_median(self, seq_gen, cutoff):
        h = smhash
        k = self.ksize
        for record in seq_gen:
            seq = dna(record.sequence)
            for kmer in sequence.kmers(self.ksize):
                counts = [min(self.counts[i, h(kmer) % p] for i, p in enumerate(self.tablesizes)) for kmer in sequence.kmers(k)]
            if np.median(counts) < cutoff:
                yield record
            else:
                count_sequence(seq)
