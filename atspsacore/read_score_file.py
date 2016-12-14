#####################################
#
#   contains method for reading a score file
#   input: file foobar.score
#   output: numpy array with arr[i][j] = score for each pair of vertices i,j
#
#####################################

import csv

import numpy as np
import pandas


def read_score_file2(file):  # deprecated, too slow
    with open(file, 'r') as f:
        f_csv = csv.reader(f, delimiter='\t')
        nr_reads = int(next(f_csv)[0])
        score = np.zeros((nr_reads, nr_reads), dtype=np.int32)
        for line in f_csv:
            i, j, val = [int(x) for x in line]
            score[i][j] = val
    return nr_reads, score


def read_score_file(file):
    with open(file, 'r') as f:
        nr_reads = int(f.readline())
    df = pandas.read_csv(file, dtype=int, skiprows=1, sep="\t", header=None, names=['i', 'j', 'score'])
    score = np.zeros((nr_reads, nr_reads), dtype=np.int32)
    for row in df.values:
        score[row[0]][row[1]] = row[2]
    return nr_reads, score

    # file = '/home/andreas/GDrive/workspace/sparsedata/ref2_c20_l100/calign.score'
    # foo1 = read_score_file(file)[1]
    # foo2 = read_score_file2(file)[1]
    # print(np.array_equal(foo1, foo2))
