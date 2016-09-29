#####################################
#
#   contains method for reading a score file
#   input: file foobar.score
#   output: numpy array with arr[i][j] = score for each pair of vertices i,j
#
#####################################

import csv

import numpy as np


def read_score_file(file):
    with open(file, 'r') as f:
        f_csv = csv.reader(f, delimiter='\t')
        nr_reads = int(next(f_csv)[0])
        score = np.zeros((nr_reads, nr_reads), dtype=np.int32)
        for line in f_csv:
            # i, j, val = [int(elem) for elem in line.split('\t')]
            i, j, val = [int(x) for x in line]
            score[i][j] = val
    return nr_reads, score
