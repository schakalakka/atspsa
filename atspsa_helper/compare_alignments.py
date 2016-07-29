import numpy as np

from config import *


def read_file(file, foo, k):
    with open(file, 'r') as f:
        f.readline()
        for line in f.readlines():
            i, j, val = [int(elem) for elem in line.split('\t')]
            foo[i][j][k] = val


def read_files(*args):
    with open(args[0], 'r') as f:
        nr_reads = int(f.readline())
    foo = np.zeros((nr_reads, nr_reads, len(args)), dtype=np.uint16)
    for k, arg in enumerate(args):
        read_file(arg, foo, k)
    return foo


for ref in [2]:
    for coverage in [20]:
        for length in [100]:
            # foo1 = read_file(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}_banded5.seqscore'.format(ref, coverage, length))
            # foo2 = read_file(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.calignscore'.format(ref, coverage, length))
            read_files(
                DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}_banded5.seqscore'.format(ref, coverage, length),
                DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.calignscore'.format(ref, coverage, length))
