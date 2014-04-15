#!/usr/bin/env python

from dna_sequence import DNASequence as dna
from itertools import izip
from smhasher import murmur3_x86_64 as smhash

class NucleotideException(Exception):
    def __init__(self, msg):
        self.msg = msg
    
    def __str__():
        return 'NucleotideError: ' + self.msg

def base_to_bits(base):
    if base == 'A':
        return 0L
    elif base == 'T':
        return 1L
    elif base == 'G':
        return 2L
    elif base == 'C':
        return 3L
    else:
        raise NucleotideException('{} not a nucleotide'.format(base))

def bits_to_base(n):
    if n == 0L:
        return 'A'
    elif n == 1L:
        return 'T'
    elif n == 2L:
        return 'G'
    elif n == 3L:
        return 'C'
    else:
        raise NucleotideException('No nucleotide for value {}'.format(n))

def comp(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    else:
        raise NucleotideException('{} is not a nucleotide!'.format(base))

def smhash_forward(kmer):

    return smhash(kmer), 0

def hash_forward(kmer):
    h = 0L
    r = 0L

    hbits = kmer.bits()
    rbits = kmer.revcomp_bits()

    h = h | hbits.next()
    r = r | rbits.next()

    for hb, rb in izip(hbits, rbits):
        h = h << 2L
        r = r << 2L
        h = h | hb
        r = r | rb

    return h, r

def hash_reverse(hashval, K):
    val = long(hashval) & 3L
    bases = [''] * K
    bases[-1] = bits_to_base(val)
    hashval = long(hashval)

    for i in xrange(1, K):
        hashval = hashval >> 2L
        val = hashval & 3L
        bases[-i-1] = bits_to_base(val)

    return dna(''.join(bases))

def main():
    print 'testing hash functions...'
    kmer1 = dna('ATCGATCGATCGATCG')
    kmer2 = dna('ATCGATCGATCGATCG')
    K = 16
    assert hash_forward(kmer1) == hash_forward(kmer2)

    h1, r1 = hash_forward(kmer1)
    h2, r2 = hash_forward(kmer2)
    assert kmer1 == hash_reverse(h1, K)
    assert hash_reverse(r1, K) == 'CGATCGATCGATCGAT'
    assert hash_reverse(r1, K) == hash_reverse(r2, K)
    assert hash_reverse(h1, K) == hash_reverse(h2, K)

if __name__ == '__main__':
    main()
