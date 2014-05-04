#!/usr/bin/env python

from smhasher import murmur3_x64_128 as smhash

class DNASequence(str):
    
    def bits(self):
        for base in self:
            if base == 'A':
                yield 0L
            elif base == 'T':
                yield 1L
            elif base == 'G':
                yield 2L
            else:
                yield 3L

    def revcomp(self):
        for base in reversed(self):
            if base == 'A':
                yield 'T'
            elif base == 'T':
                yield 'A'
            elif base == 'G':
                yield 'C'
            else:
                yield 'G'

    def revcomp_bits(self):
        for base in reversed(self):
            if base == 'A':
                yield 1L
            elif base == 'T':
                yield 0L
            elif base == 'G':
                yield 3L
            else:
                yield 2L

    def kmers(self, K):
        for i in xrange(len(self)-K+1):
            yield DNASequence(self[i:i+K])

    def kmers_hashed(self, K):
        h = smhash
        for i in xrange(len(self)-K+1):
            yield h(self[i:i+K])
