#!/usr/bin/env python

import array
import sys

from dna_sequence import DNASequence as dna
from pykhmer_utils import get_n_primes_above_x as get_primes
from pykhmer_utils import KMerException, NucleotideException, zeros
from smhasher import murmur3_x64_128 as smhash
import numpy as np


class KCMSException(Exception):
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return 'KMerCountMinSketchError: ' + self.msg

class KMerCountMinSketch(object):

    def __init__(self, ksize, minhashsize, nhashes, dtype=None, arrtype=bytearray):
        
        self.valid_types = [int, np.uint64, np.uint8, np.uint16, np.uint32, None]
        
        if ksize < 5:
            raise KCMSException('ksize too small:' + ksize)
        if nhashes < 1:
            raise KCMSException('nhashes must be > 0')
        if minhashsize < 1e5:
            raise KCMSException('minhashsize should be > 10000')
        
        self.ksize = ksize
        self.minhashsize = long(minhashsize)
        self.nhashes = nhashes

        self.dtype = dtype
        self.arrtype = arrtype
        
        self.primes = get_primes(self.nhashes, self.minhashsize)
        
        if arrtype is bytearray:
            if dtype != None:
                raise TypeError('Cannot specific dtype when arrtype is bytearray')
            self.counts = [bytearray(zeros(p)) for p in self.primes]
        else:
            if dtype not in self.valid_types:
                raise TypeError('dtype must be ' + ', '.join(str(t) for t in self.valid_types))
            if arrtype is array.array:
                if self.dtype in [int, np.uint64, np.uint32]:
                    self.counts = [array.array('L', zeros(p, dtype=self.dtype)) for p in self.primes]
                elif self.dtype in [np.uint16]:
                    self.counts = [array.array('H', zeros(p, dtype=self.dtype)) for p in self.primes]
                else:
                    self.counts = [array.array('B', zeros(p)) for p in self.primes]
            elif arrtype is np.array:
                if self.dtype == None:
                    self.counts = [np.array([z for z in zeros(p)]) for p in self.primes]
                else:
                    self.counts = [np.array([z for z in zeros(p)], dtype=self.dtype) for p in self.primes]
            else:
                # this could definitely cause exceptions
                self.counts = [arrtype([z for z in zeros(p)]) for p in self.primes]
        
        for table in self.counts:
            assert len(table) >= 1e5
            for c in table:
                assert c == 0
        assert self.ksize > 4
        assert len(self.counts) == self.nhashes
        assert set(self.primes) == set([len(t) for t in self.counts])
        
        print >>sys.stderr, 'KMerCountMinSketch\n'+'*'*18
        print >>sys.stderr, 'K:', self.ksize
        print >>sys.stderr, 'nhashes:', self.nhashes
        print >>sys.stderr, 'Hashing primes:\n', ', '.join(str(t) for t in self.primes)
        print >>sys.stderr, 'Array type:', type(self.counts[0])
        print >>sys.stderr, 'Array dtype:', type(self.counts[0][0])

        self.median = np.median

    def get(self, kmer):
        if len(kmer) != self.ksize:
            raise KMerException('Tried to get count \
                                with invalid k-mer size')

        kh = smhash(kmer)
        counts = self.counts
        return min(table[kh % len(table)] for table in counts)
    
    def get_hashed(self, kh):
        counts = self.counts
        return min(table[kh % len(table)] for table in counts)

    def count(self, kmer):
        kh = smhash(kmer)
        counts = self.counts
        for i, p in enumerate(self.primes):
            try:
                counts[i][kh % p] += 1
            except (ValueError, OverflowError):
                pass

    def count_sequence(self, sequence):
        primes = self.primes
        counts = self.counts
        for kh in sequence.kmers_hashed(self.ksize):
            for i, p in enumerate(primes):
                try:
                    counts[i][kh % p] += 1
                except (ValueError, OverflowError):
                    pass

    def n_occupied(self):
        return min(np.count_nonzero(row) for row in self.counts)

    def get_median_count(self, sequence):
        get = self.get_hashed
        return self.median([get(kh) for kh in sequence.kmers_hashed(self.ksize)])
    
    # sequences should be a generator that produces (name, dna(sequence)) tuples
    def normalize_by_median(self, sequences, cutoff):
        k = self.ksize
        primes = self.primes
        counts = self.counts
        for name, sequence in sequences:
            kmer_counts = [min(counts[i][kh % p] for i, p in enumerate(primes)) for kh in sequence.kmers_hashed(k)]
            if np.median(kmer_counts) < cutoff:
                yield name, sequence
                for kh in sequence.kmers_hashed(self.ksize):
                    for i, p in enumerate(primes):
                        try:
                            counts[i][kh % p] += 1
                        except (ValueError, OverflowError):
                            pass
