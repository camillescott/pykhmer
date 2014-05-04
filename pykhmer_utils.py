#!/usr/bin/env python

from dna_sequence import DNASequence as dna
import random

bases = ['A', 'T', 'C', 'G']
def rand_dna(L):
    s = bases * (L // 4)
    random.shuffle(s)
    return dna(''.join(s))

def rand_dnas(n, L=100):
    r = rand_dna
    for i in xrange(n):
        yield rand_dna(L)

def fasta_generator(n, L):
    for i, seq in enumerate(rand_dnas(n, L)):
        yield '>Seq {}'.format(i), seq

class KMerException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return 'KMerException: ' + self.msg

class NucleotideException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__():
        return 'NucleotideError: ' + self.msg

def zeros(n, dtype=np.uint8):
    while(n > 0):
        yield dtype(0)
        n = n - 1

# utilities, some borrowed from khmer

def calc_expected_collisions(hashtable):
    """Do a quick & dirty expected collision rate calculation on a hashtable.

    Keyword argument:
    hashtable: the hashtable object to inspect
    """
    sizes = hashtable.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(hashtable.n_occupied())
    min_size = min(sizes)

    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht

    return fp_all


def is_prime(number):
    '''Checks if a number is prime.'''
    if number < 2:
        return False
    if number == 2:
        return True
    if number % 2 == 0:
        return False
    for _ in range(3, int(number ** 0.5) + 1, 2):
        if number % _ == 0:
            return False
    return True


def get_n_primes_near_x(number, target):
    ''' Step backwards until a number of primes (other than 2) have been
    found that are smaller than the target and return them.

    Keyword arguments:
    number -- the number of primes to find
    target -- the number to step backwards from
    '''
    primes = []
    i = target - 1
    if i % 2 == 0:
        i -= 1
    while len(primes) != number and i > 0:
        if is_prime(i):
            primes.append(i)
        i -= 2
    return primes


def get_n_primes_above_x(number, target):
    '''Step forwards until a number of primes (other than 2) have been
    found that are smaller than the target and return them.

    Keyword arguments:
    number -- the number of primes to find
    target -- the number to step forwards from
    '''
    primes = []
    i = target + 1
    if i % 2 == 0:
        i += 1
    while len(primes) != number and i > 0:
        if is_prime(i):
            primes.append(i)
        i += 2
    return primes


