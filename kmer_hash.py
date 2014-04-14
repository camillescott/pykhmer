#!/usr/bin/env python

class NucleotideException(Exception):
    def __init__(self, msg):
        self.msg = msg
    
    def __str__():
        return 'NucleotideError: ' + self.msg

def base_to_bits(base):
    if base == 'A' or base == 'a':
        return 0L
    elif base == 'T' or base == 't':
        return 1L
    elif base == 'G' or base == 'g':
        return 2L
    elif base == 'C' or base == 'c':
        return 3L
    else:
        raise NucleotideException('{} is not a nucleotide!'.format(base))

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

def hash_forward(kmer, K):
    h = 0L
    r = 0L

    h = h | base_to_bits(kmer[0])
    r = r | base_to_bits(comp(kmer[-1]))

    for i in xrange(1, K):
        j = K - 1 - i
        h = h << 2L
        r = r << 2L
        h = h | base_to_bits(kmer[i])
        r = r | base_to_bits(comp(kmer[j]))
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

    return ''.join(bases)

def main():
    print 'testing hash functions...'
    kmer1 = 'ATCGATCGATCGATCG'
    kmer2 = 'ATCGATCGATCGATCG'
    K = 16
    assert hash_forward(kmer1, K) == hash_forward(kmer2, K)

    h1, r1 = hash_forward(kmer1, K)
    h2, r2 = hash_forward(kmer2, K)
    assert kmer1 == hash_reverse(h1, K)
    assert hash_reverse(r1, K) == 'CGATCGATCGATCGAT'
    assert hash_reverse(r1, K) == hash_reverse(r2, K)
    assert hash_reverse(h1, K) == hash_reverse(h2, K)

if __name__ == '__main__':
    main()
